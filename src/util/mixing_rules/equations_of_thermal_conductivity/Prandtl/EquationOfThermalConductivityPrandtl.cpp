#include "util/mixing_rules/equations_of_thermal_conductivity/Prandtl/EquationOfThermalConductivityPrandtl.hpp"

/*
 * Print all characteristics of the equation of thermal conductivity class.
 */
void
EquationOfThermalConductivityPrandtl::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfThermalConductivityPrandtl object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfThermalConductivityPrandtl: this = "
       << (EquationOfThermalConductivityPrandtl *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the thermal conductivity.
 */
double
EquationOfThermalConductivityPrandtl::getThermalConductivity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& molecular_properties) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 4);
#endif
    
    const double& c_p = *(molecular_properties[0]);
    const double& Pr = *(molecular_properties[1]);
    
    std::vector<double> mu_molecular_properties;
    std::vector<const double*> mu_molecular_properties_const_ptr;
    
    mu_molecular_properties.reserve(static_cast<int>(molecular_properties.size()) - 3);
    mu_molecular_properties_const_ptr.reserve(static_cast<int>(molecular_properties.size()) - 3);
    
    for (int mi = 3; mi < static_cast<int>(molecular_properties.size()); mi++)
    {
        mu_molecular_properties.push_back(*molecular_properties[mi]);
    }
    
    for (int mi = 0; mi < (static_cast<int>(molecular_properties.size()) - 3); mi++)
    {
        mu_molecular_properties_const_ptr.push_back(&mu_molecular_properties[mi]);
    }
    
    const double mu = d_equation_of_shear_viscosity->
        getShearViscosity(
            pressure,
            temperature,
            mu_molecular_properties_const_ptr);
    
    return c_p*mu/Pr;
}


/*
 * Compute the thermal conductivity.
 */
void
EquationOfThermalConductivityPrandtl::computeThermalConductivity(
    boost::shared_ptr<pdat::CellData<double> >& data_thermal_conductivity,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const std::vector<const double*>& molecular_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_thermal_conductivity);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 4);
#endif
    
    // Get the dimensions of the ghost cell box.
    const hier::Box ghost_box_thermal_conductivity = data_thermal_conductivity->getGhostBox();
    const hier::IntVector ghostcell_dims_thermal_conductivity = ghost_box_thermal_conductivity.numberCells();
    
    // Delcare data container for shear viscosity.
    boost::shared_ptr<pdat::CellData<double> > data_shear_viscosity;
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets of all data and dimensions of the ghost cell box for shear viscosity
     * and allocate memory.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_thermal_conductivity(d_dim);
    hier::IntVector offset_min(d_dim);
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_thermal_conductivity = data_thermal_conductivity->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        
        // Get the interior box and the dimensions of box that covers the interior of patch.
        const hier::Box interior_box = data_thermal_conductivity->getBox();
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
        
        num_ghosts_min = num_ghosts_thermal_conductivity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_min = num_ghosts_min;
        offset_thermal_conductivity = num_ghosts_thermal_conductivity;
        
        ghostcell_dims_min = interior_dims + num_ghosts_min*2;
        
        data_shear_viscosity = boost::make_shared<pdat::CellData<double> >(interior_box, 1, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_thermal_conductivity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_min = hier::IntVector::getZero(d_dim);
        offset_thermal_conductivity = domain.lower() - ghost_box_thermal_conductivity.lower();
        
        ghostcell_dims_min = domain_dims;
        
        data_shear_viscosity = boost::make_shared<pdat::CellData<double> >(domain, 1, hier::IntVector::getZero(d_dim));
    }
    
    /*
     * Compute the shear viscosity.
     */
    
    std::vector<double> mu_molecular_properties;
    std::vector<const double*> mu_molecular_properties_const_ptr;
    
    mu_molecular_properties.reserve(static_cast<int>(molecular_properties.size()) - 3);
    mu_molecular_properties_const_ptr.reserve(static_cast<int>(molecular_properties.size()) - 3);
    
    for (int mi = 3; mi < static_cast<int>(molecular_properties.size()); mi++)
    {
        mu_molecular_properties.push_back(*molecular_properties[mi]);
    }
    
    for (int mi = 0; mi < (static_cast<int>(molecular_properties.size()) - 3); mi++)
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
    
    double* kappa = data_thermal_conductivity->getPointer(0);
    double* mu = data_shear_viscosity->getPointer(0);
    
    const double& c_p = *(molecular_properties[0]);
    const double& Pr = *(molecular_properties[1]);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_thermal_conductivity = offset_thermal_conductivity[0];
        const int offset_0_min = offset_min[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_thermal_conductivity = i + offset_0_thermal_conductivity;
            const int idx_min = i + offset_0_min;
            
            kappa[idx_thermal_conductivity] = c_p*mu[idx_min]/Pr;
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
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_thermal_conductivity = (i + offset_0_thermal_conductivity) +
                    (j + offset_1_thermal_conductivity)*ghostcell_dim_0_thermal_conductivity;
                
                const int idx_min = (i + offset_0_min) +
                    (j + offset_1_min)*ghostcell_dim_0_min;
                
                kappa[idx_thermal_conductivity] = c_p*mu[idx_min]/Pr;
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
                    const int idx_thermal_conductivity = (i + offset_0_thermal_conductivity) +
                        (j + offset_1_thermal_conductivity)*ghostcell_dim_0_thermal_conductivity +
                        (k + offset_2_thermal_conductivity)*ghostcell_dim_0_thermal_conductivity*
                            ghostcell_dim_1_thermal_conductivity;
                    
                    const int idx_min = (i + offset_0_min) +
                        (j + offset_1_min)*ghostcell_dim_0_min +
                        (k + offset_2_min)*ghostcell_dim_0_min*
                            ghostcell_dim_1_min;
                    
                    kappa[idx_thermal_conductivity] = c_p*mu[idx_min]/Pr;
                }
            }
        }
    }
}


/*
 * Compute the thermal conductivity.
 */
void
EquationOfThermalConductivityPrandtl::computeThermalConductivity(
    boost::shared_ptr<pdat::CellData<double> >& data_thermal_conductivity,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_molecular_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_thermal_conductivity);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_molecular_properties);
    
    TBOX_ASSERT(data_molecular_properties->getDepth() >= 4);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_thermal_conductivity = data_thermal_conductivity->getGhostBox();
    const hier::IntVector ghostcell_dims_thermal_conductivity = ghost_box_thermal_conductivity.numberCells();
    
    const hier::Box ghost_box_molecular_properties = data_molecular_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_molecular_properties = ghost_box_molecular_properties.numberCells();
    
    // Delcare data containers for shear viscosity and molecular properties.
    boost::shared_ptr<pdat::CellData<double> > data_shear_viscosity;
    boost::shared_ptr<pdat::CellData<double> > data_molecular_properties_shear_viscosity;
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets of all data and dimensions of the ghost cell box for shear viscosity
     * and molecular properties and allocate memory.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_thermal_conductivity(d_dim);
    hier::IntVector offset_molecular_properties(d_dim);
    hier::IntVector offset_min(d_dim);
    
    hier::IntVector ghostcell_dims_min(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_thermal_conductivity = data_thermal_conductivity->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_molecular_properties = data_molecular_properties->getGhostCellWidth();
        
        // Get the interior box and the dimensions of box that covers the interior of patch.
        const hier::Box interior_box = data_thermal_conductivity->getBox();
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
        
        num_ghosts_min = num_ghosts_thermal_conductivity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_molecular_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_min = num_ghosts_min;
        offset_thermal_conductivity = num_ghosts_thermal_conductivity;
        offset_molecular_properties = num_ghosts_molecular_properties;
        
        ghostcell_dims_min = interior_dims + num_ghosts_min*2;
        
        data_shear_viscosity = boost::make_shared<pdat::CellData<double> >(interior_box, 1, num_ghosts_min);
        
        data_molecular_properties_shear_viscosity = boost::make_shared<pdat::CellData<double> >(
            interior_box, data_molecular_properties->getDepth() - 3, num_ghosts_min);
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_thermal_conductivity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_molecular_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_min = hier::IntVector::getZero(d_dim);
        offset_thermal_conductivity = domain.lower() - ghost_box_thermal_conductivity.lower();
        offset_molecular_properties = domain.lower() - ghost_box_molecular_properties.lower();
        
        ghostcell_dims_min = domain_dims;
        
        data_shear_viscosity = boost::make_shared<pdat::CellData<double> >(domain, 1, hier::IntVector::getZero(d_dim));
        
        data_molecular_properties_shear_viscosity = boost::make_shared<pdat::CellData<double> >(
            domain, data_molecular_properties->getDepth() - 3, hier::IntVector::getZero(d_dim));
    }
    
    /*
     * Compute the shear viscosity.
     */
    
    for (int mi = 0; mi < data_molecular_properties->getDepth() - 3; mi++)
    {
        data_molecular_properties_shear_viscosity->copyDepth(mi, *data_molecular_properties, mi + 3);
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
    
    double* kappa = data_thermal_conductivity->getPointer(0);
    double* mu = data_shear_viscosity->getPointer(0);
    
    double* c_p = data_molecular_properties->getPointer(0);
    double* Pr = data_molecular_properties->getPointer(1);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_thermal_conductivity = offset_thermal_conductivity[0];
        const int offset_0_molecular_properties = offset_molecular_properties[0];
        const int offset_0_min = offset_min[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_thermal_conductivity = i + offset_0_thermal_conductivity;
            const int idx_molecular_properties = i + offset_0_molecular_properties;
            const int idx_min = i + offset_0_min;
            
            kappa[idx_thermal_conductivity] = c_p[idx_molecular_properties]*mu[idx_min]/
                Pr[idx_molecular_properties];
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
        
        const int offset_0_molecular_properties = offset_molecular_properties[0];
        const int offset_1_molecular_properties = offset_molecular_properties[1];
        const int ghostcell_dim_0_molecular_properties = ghostcell_dims_molecular_properties[0];
        
        const int offset_0_min = offset_min[0];
        const int offset_1_min = offset_min[1];
        const int ghostcell_dim_0_min = ghostcell_dims_min[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_thermal_conductivity = (i + offset_0_thermal_conductivity) +
                    (j + offset_1_thermal_conductivity)*ghostcell_dim_0_thermal_conductivity;
                
                const int idx_molecular_properties = (i + offset_0_molecular_properties) +
                    (j + offset_1_molecular_properties)*ghostcell_dim_0_molecular_properties;
                
                const int idx_min = (i + offset_0_min) +
                    (j + offset_1_min)*ghostcell_dim_0_min;
                
                kappa[idx_thermal_conductivity] = c_p[idx_molecular_properties]*mu[idx_min]/
                    Pr[idx_molecular_properties];
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
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_thermal_conductivity = (i + offset_0_thermal_conductivity) +
                        (j + offset_1_thermal_conductivity)*ghostcell_dim_0_thermal_conductivity +
                        (k + offset_2_thermal_conductivity)*ghostcell_dim_0_thermal_conductivity*
                            ghostcell_dim_1_thermal_conductivity;
                    
                    const int idx_molecular_properties = (i + offset_0_molecular_properties) +
                        (j + offset_1_molecular_properties)*ghostcell_dim_0_molecular_properties +
                        (k + offset_2_molecular_properties)*ghostcell_dim_0_molecular_properties*
                            ghostcell_dim_1_molecular_properties;
                    
                    const int idx_min = (i + offset_0_min) +
                        (j + offset_1_min)*ghostcell_dim_0_min +
                        (k + offset_2_min)*ghostcell_dim_0_min*
                            ghostcell_dim_1_min;
                    
                    kappa[idx_thermal_conductivity] = c_p[idx_molecular_properties]*mu[idx_min]/
                        Pr[idx_molecular_properties];
                }
            }
        }
    }
}
