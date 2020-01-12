#include "util/mixing_rules/equations_of_shear_viscosity/constant/EquationOfShearViscosityConstant.hpp"

/*
 * Print all characteristics of the equation of shear viscosity class.
 */
void
EquationOfShearViscosityConstant::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfShearViscosityConstant object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfShearViscosityConstant: this = "
       << (EquationOfShearViscosityConstant *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the shear viscosity.
 */
double
EquationOfShearViscosityConstant::getShearViscosity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& molecular_properties) const
{
    NULL_USE(pressure);
    NULL_USE(temperature);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 1);
#endif
    
    return *(molecular_properties[0]);
}


/*
 * Compute the shear viscosity.
 */
void
EquationOfShearViscosityConstant::computeShearViscosity(
    boost::shared_ptr<pdat::CellData<double> >& data_shear_viscosity,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const std::vector<const double*>& molecular_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_pressure);
    NULL_USE(data_temperature);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_shear_viscosity);
    
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 1);
#endif
    
    if (domain.empty())
    {
        data_shear_viscosity->fillAll(*molecular_properties[0]);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_shear_viscosity->getGhostBox().contains(domain));
#endif
        
        data_shear_viscosity->fillAll(*molecular_properties[0], domain);
    }
}


/*
 * Compute the shear viscosity.
 */
void
EquationOfShearViscosityConstant::computeShearViscosity(
    boost::shared_ptr<pdat::CellData<double> >& data_shear_viscosity,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_molecular_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_pressure);
    NULL_USE(data_temperature);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_shear_viscosity);
    TBOX_ASSERT(data_molecular_properties);
    
    TBOX_ASSERT(data_molecular_properties->getDepth() >= 1);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_shear_viscosity = data_shear_viscosity->getGhostBox();
    const hier::IntVector ghostcell_dims_shear_viscosity = ghost_box_shear_viscosity.numberCells();
    
    const hier::Box ghost_box_molecular_properties = data_molecular_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_molecular_properties = ghost_box_molecular_properties.numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_shear_viscosity(d_dim);
    hier::IntVector offset_molecular_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the number of ghost cells.
        const hier::IntVector num_ghosts_shear_viscosity = data_shear_viscosity->getGhostCellWidth();
        const hier::IntVector num_ghosts_molecular_properties = data_molecular_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_shear_viscosity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_molecular_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_shear_viscosity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_molecular_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_shear_viscosity = num_ghosts_shear_viscosity;
        offset_molecular_properties = num_ghosts_molecular_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_shear_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_molecular_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_shear_viscosity = domain.lower() - ghost_box_shear_viscosity.lower();
        offset_molecular_properties = domain.lower() - ghost_box_molecular_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* mu = data_shear_viscosity->getPointer(0);
    double* mu_src = data_molecular_properties->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_shear_viscosity = offset_shear_viscosity[0];
        const int offset_0_molecular_properties = offset_molecular_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_shear_viscosity = i + offset_0_shear_viscosity;
            const int idx_molecular_properties = i + offset_0_molecular_properties;
            
            mu[idx_shear_viscosity] = mu_src[idx_molecular_properties];
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
        
        const int offset_0_shear_viscosity = offset_shear_viscosity[0];
        const int offset_1_shear_viscosity = offset_shear_viscosity[1];
        const int ghostcell_dim_0_shear_viscosity = ghostcell_dims_shear_viscosity[0];
        
        const int offset_0_molecular_properties = offset_molecular_properties[0];
        const int offset_1_molecular_properties = offset_molecular_properties[1];
        const int ghostcell_dim_0_molecular_properties = ghostcell_dims_molecular_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_shear_viscosity = (i + offset_0_shear_viscosity) +
                    (j + offset_1_shear_viscosity)*ghostcell_dim_0_shear_viscosity;
                
                const int idx_molecular_properties = (i + offset_0_molecular_properties) +
                    (j + offset_1_molecular_properties)*ghostcell_dim_0_molecular_properties;
                
                mu[idx_shear_viscosity] = mu_src[idx_molecular_properties];
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
        
        const int offset_0_shear_viscosity = offset_shear_viscosity[0];
        const int offset_1_shear_viscosity = offset_shear_viscosity[1];
        const int offset_2_shear_viscosity = offset_shear_viscosity[2];
        const int ghostcell_dim_0_shear_viscosity = ghostcell_dims_shear_viscosity[0];
        const int ghostcell_dim_1_shear_viscosity = ghostcell_dims_shear_viscosity[1];
        
        const int offset_0_molecular_properties = offset_molecular_properties[0];
        const int offset_1_molecular_properties = offset_molecular_properties[1];
        const int offset_2_molecular_properties = offset_molecular_properties[2];
        const int ghostcell_dim_0_molecular_properties = ghostcell_dims_molecular_properties[0];
        const int ghostcell_dim_1_molecular_properties = ghostcell_dims_molecular_properties[1];
        
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
                    const int idx_shear_viscosity = (i + offset_0_shear_viscosity) +
                        (j + offset_1_shear_viscosity)*ghostcell_dim_0_shear_viscosity +
                        (k + offset_2_shear_viscosity)*ghostcell_dim_0_shear_viscosity*
                            ghostcell_dim_1_shear_viscosity;
                    
                    const int idx_molecular_properties = (i + offset_0_molecular_properties) +
                        (j + offset_1_molecular_properties)*ghostcell_dim_0_molecular_properties +
                        (k + offset_2_molecular_properties)*ghostcell_dim_0_molecular_properties*
                            ghostcell_dim_1_molecular_properties;
                    
                    mu[idx_shear_viscosity] = mu_src[idx_molecular_properties];
                }
            }
        }
    }
}
