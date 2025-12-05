#include "util/mixing_rules/equations_of_bulk_viscosity/constant_ratio_to_shear_viscosity/EquationOfBulkViscosityConstantRatioToShearViscosity.hpp"

/*
 * Print all characteristics of the equation of bulk viscosity class.
 */
void
EquationOfBulkViscosityConstantRatioToShearViscosity::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfBulkViscosityConstantRatioToShearViscosity object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfBulkViscosityConstantRatioToShearViscosity: this = "
       << (EquationOfBulkViscosityConstantRatioToShearViscosity *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the bulk viscosity.
 */
Real
EquationOfBulkViscosityConstantRatioToShearViscosity::getBulkViscosity(
    const Real* const pressure,
    const Real* const temperature,
    const std::vector<const Real*>& molecular_properties) const
{
    NULL_USE(pressure);
    NULL_USE(temperature);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 3);
#endif
    
    const Real& ratio_of_bulk_viscosity_to_shear_viscosity = *(molecular_properties[0]);
    
    std::vector<Real> mu_molecular_properties;
    std::vector<const Real*> mu_molecular_properties_const_ptr;
    
    mu_molecular_properties.reserve(static_cast<int>(molecular_properties.size()) - 2);
    mu_molecular_properties_const_ptr.reserve(static_cast<int>(molecular_properties.size()) - 2);
    
    for (int mi = 2; mi < static_cast<int>(molecular_properties.size()); mi++)
    {
        mu_molecular_properties.push_back(*molecular_properties[mi]);
    }
    
    for (int mi = 0; mi < (static_cast<int>(molecular_properties.size()) - 2); mi++)
    {
        mu_molecular_properties_const_ptr.push_back(&mu_molecular_properties[mi]);
    }
    
    const Real mu = d_equation_of_shear_viscosity->
        getShearViscosity(
            pressure,
            temperature,
            mu_molecular_properties_const_ptr);
    
    return ratio_of_bulk_viscosity_to_shear_viscosity*mu;
}


/*
 * Compute the bulk viscosity.
 */
void
EquationOfBulkViscosityConstantRatioToShearViscosity::computeBulkViscosity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_bulk_viscosity,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const std::vector<const Real*>& molecular_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_bulk_viscosity);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 3);
#endif
    
    // Get the dimensions of the ghost cell box.
    const hier::Box ghost_box_bulk_viscosity = data_bulk_viscosity->getGhostBox();
    const hier::IntVector ghostcell_dims_bulk_viscosity = ghost_box_bulk_viscosity.numberCells();
    
    // Delcare data container for shear viscosity.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_shear_viscosity;
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets of all data and dimensions of the ghost cell box for shear viscosity
     * and allocate memory.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_bulk_viscosity(d_dim);
    hier::IntVector offset_min(d_dim);
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_bulk_viscosity = data_bulk_viscosity->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        
        // Get the interior box and the dimensions of box that covers the interior of patch.
        const hier::Box interior_box = data_bulk_viscosity->getBox();
        const hier::IntVector interior_dims = interior_box.numberCells();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells and the dimensions of the ghost cell box for shear
         * viscosity.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_bulk_viscosity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_min = num_ghosts_min;
        offset_bulk_viscosity = num_ghosts_bulk_viscosity;
        
        ghostcell_dims_min = interior_dims + num_ghosts_min*2;
        
        data_shear_viscosity = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(interior_box, 1, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_bulk_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_min = hier::IntVector::getZero(d_dim);
        offset_bulk_viscosity = domain.lower() - ghost_box_bulk_viscosity.lower();
        
        ghostcell_dims_min = domain_dims;
        
        data_shear_viscosity = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(domain, 1, hier::IntVector::getZero(d_dim));
    }
    
    /*
     * Compute the shear viscosity.
     */
    
    std::vector<Real> mu_molecular_properties;
    std::vector<const Real*> mu_molecular_properties_const_ptr;
    
    mu_molecular_properties.reserve(static_cast<int>(molecular_properties.size()) - 2);
    mu_molecular_properties_const_ptr.reserve(static_cast<int>(molecular_properties.size()) - 2);
    
    for (int mi = 2; mi < static_cast<int>(molecular_properties.size()); mi++)
    {
        mu_molecular_properties.push_back(*molecular_properties[mi]);
    }
    
    for (int mi = 0; mi < (static_cast<int>(molecular_properties.size()) - 2); mi++)
    {
        mu_molecular_properties_const_ptr.push_back(&mu_molecular_properties[mi]);
    }
    
    d_equation_of_shear_viscosity->computeShearViscosity(
        data_shear_viscosity,
        data_pressure,
        data_temperature,
        mu_molecular_properties_const_ptr,
        domain);
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* mu_v = data_bulk_viscosity->getPointer(0);
    Real* mu = data_shear_viscosity->getPointer(0);
    
    const Real& ratio_of_bulk_viscosity_to_shear_viscosity = *(molecular_properties[0]);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_bulk_viscosity = offset_bulk_viscosity[0];
        const int offset_0_min = offset_min[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_bulk_viscosity = i + offset_0_bulk_viscosity;
            const int idx_min = i + offset_0_min;
            
            mu_v[idx_bulk_viscosity] = ratio_of_bulk_viscosity_to_shear_viscosity*mu[idx_min];
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
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                    (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity;
                
                const int idx_min = (i + offset_0_min) +
                    (j + offset_1_min)*ghostcell_dim_0_min;
                
                mu_v[idx_bulk_viscosity] = ratio_of_bulk_viscosity_to_shear_viscosity*mu[idx_min];
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
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    mu_v[idx_bulk_viscosity] = ratio_of_bulk_viscosity_to_shear_viscosity*mu[idx_min];
                }
            }
        }
    }
}


/*
 * Compute the bulk viscosity.
 */
void
EquationOfBulkViscosityConstantRatioToShearViscosity::computeBulkViscosity(
    HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_bulk_viscosity,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_molecular_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_bulk_viscosity);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_molecular_properties);
    
    TBOX_ASSERT(data_molecular_properties->getDepth() >= 3);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_bulk_viscosity = data_bulk_viscosity->getGhostBox();
    const hier::IntVector ghostcell_dims_bulk_viscosity = ghost_box_bulk_viscosity.numberCells();
    
    const hier::Box ghost_box_molecular_properties = data_molecular_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_molecular_properties = ghost_box_molecular_properties.numberCells();
    
    // Delcare data containers for shear viscosity and molecular properties.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_shear_viscosity;
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_molecular_properties_shear_viscosity;
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets of all data and dimensions of the ghost cell box for shear viscosity
     * and molecular properties and allocate memory.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_bulk_viscosity(d_dim);
    hier::IntVector offset_molecular_properties(d_dim);
    hier::IntVector offset_min(d_dim);
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_bulk_viscosity = data_bulk_viscosity->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_molecular_properties = data_molecular_properties->getGhostCellWidth();
        
        // Get the interior box and the dimensions of box that covers the interior of patch.
        const hier::Box interior_box = data_bulk_viscosity->getBox();
        const hier::IntVector interior_dims = interior_box.numberCells();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_molecular_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        /*
         * Get the minimum number of ghost cells and the dimensions of the ghost cell box for shear
         * viscosity and molecular properties.
         */
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_bulk_viscosity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_molecular_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_min = num_ghosts_min;
        offset_bulk_viscosity = num_ghosts_bulk_viscosity;
        offset_molecular_properties = num_ghosts_molecular_properties;
        
        ghostcell_dims_min = interior_dims + num_ghosts_min*2;
        
        data_shear_viscosity = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(interior_box, 1, num_ghosts_min);
        
        data_molecular_properties_shear_viscosity = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
            interior_box, data_molecular_properties->getDepth() - 2, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_bulk_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_molecular_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_min = hier::IntVector::getZero(d_dim);
        offset_bulk_viscosity = domain.lower() - ghost_box_bulk_viscosity.lower();
        offset_molecular_properties = domain.lower() - ghost_box_molecular_properties.lower();
        
        ghostcell_dims_min = domain_dims;
        
        data_shear_viscosity = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(domain, 1, hier::IntVector::getZero(d_dim));
        
        data_molecular_properties_shear_viscosity = HAMERS_MAKE_SHARED<pdat::CellData<Real> >(
            domain, data_molecular_properties->getDepth() - 2, hier::IntVector::getZero(d_dim));
    }
    
    /*
     * Compute the shear viscosity.
     */
    
    for (int mi = 0; mi < data_molecular_properties->getDepth() - 2; mi++)
    {
        data_molecular_properties_shear_viscosity->copyDepth(mi, *data_molecular_properties, mi + 2);
    }
    
    d_equation_of_shear_viscosity->computeShearViscosity(
        data_shear_viscosity,
        data_pressure,
        data_temperature,
        data_molecular_properties_shear_viscosity,
        domain);
    
    /*
     * Get the pointers to the cell data.
     */
    
    Real* mu_v = data_bulk_viscosity->getPointer(0);
    Real* mu = data_shear_viscosity->getPointer(0);
    
    Real* ratio_of_bulk_viscosity_to_shear_viscosity = data_molecular_properties->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_bulk_viscosity = offset_bulk_viscosity[0];
        const int offset_0_molecular_properties = offset_molecular_properties[0];
        const int offset_0_min = offset_min[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_bulk_viscosity = i + offset_0_bulk_viscosity;
            const int idx_molecular_properties = i + offset_0_molecular_properties;
            const int idx_min = i + offset_0_min;
            
            mu_v[idx_bulk_viscosity] = ratio_of_bulk_viscosity_to_shear_viscosity[idx_molecular_properties]*mu[idx_min];
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
        
        const int offset_0_molecular_properties = offset_molecular_properties[0];
        const int offset_1_molecular_properties = offset_molecular_properties[1];
        const int ghostcell_dim_0_molecular_properties = ghostcell_dims_molecular_properties[0];
        
        const int offset_0_min = offset_min[0];
        const int offset_1_min = offset_min[1];
        const int ghostcell_dim_0_min = ghostcell_dims_min[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                    (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity;
                
                const int idx_molecular_properties = (i + offset_0_molecular_properties) +
                    (j + offset_1_molecular_properties)*ghostcell_dim_0_molecular_properties;
                
                const int idx_min = (i + offset_0_min) +
                    (j + offset_1_min)*ghostcell_dim_0_min;
                
                mu_v[idx_bulk_viscosity] = ratio_of_bulk_viscosity_to_shear_viscosity[idx_molecular_properties]*mu[idx_min];
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
        
        const int offset_0_molecular_properties = offset_molecular_properties[0];
        const int offset_1_molecular_properties = offset_molecular_properties[1];
        const int offset_2_molecular_properties = offset_molecular_properties[2];
        const int ghostcell_dim_0_molecular_properties = ghostcell_dims_molecular_properties[0];
        const int ghostcell_dim_1_molecular_properties = ghostcell_dims_molecular_properties[1];
        
        const int offset_0_min = offset_min[0];
        const int offset_1_min = offset_min[1];
        const int offset_2_min = offset_min[2];
        const int ghostcell_dim_0_min = ghostcell_dims_min[0];
        const int ghostcell_dim_1_min = ghostcell_dims_min[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                        (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity +
                        (k + offset_2_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity*
                            ghostcell_dim_1_bulk_viscosity;
                    
                    const int idx_molecular_properties = (i + offset_0_molecular_properties) +
                        (j + offset_1_molecular_properties)*ghostcell_dim_0_molecular_properties +
                        (k + offset_2_molecular_properties)*ghostcell_dim_0_molecular_properties*
                            ghostcell_dim_1_molecular_properties;
                    
                    const int idx_min = (i + offset_0_min) +
                        (j + offset_1_min)*ghostcell_dim_0_min +
                        (k + offset_2_min)*ghostcell_dim_0_min*
                            ghostcell_dim_1_min;
                    
                    mu_v[idx_bulk_viscosity] = ratio_of_bulk_viscosity_to_shear_viscosity[idx_molecular_properties]*mu[idx_min];
                }
            }
        }
    }
}
