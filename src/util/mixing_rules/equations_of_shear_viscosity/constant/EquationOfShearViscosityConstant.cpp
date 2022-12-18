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
    // If the constant kinematic viscosity is used.
    if (d_use_constant_kinematic_viscosity_and_ideal_gas_assumptions)
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 2);
#endif
        
        const double& nu = *(molecular_properties[0]);
        const double& M  = *(molecular_properties[1]);
        
        const double R = d_R_u/M;
        
        const double& p = *pressure;
        const double& T = *temperature;
        
        const double rho = p/(R*T);
        
        return (rho*nu);
    }
    // If the constant kinematic dynamic is used.
    else
    {
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
        TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 1);
#endif
        
        return *(molecular_properties[0]);
    }
}


/*
 * Compute the shear viscosity.
 */
void
EquationOfShearViscosityConstant::computeShearViscosity(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_shear_viscosity,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
    const std::vector<const double*>& molecular_properties,
    const hier::Box& domain) const
{

#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_shear_viscosity);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    
    // If the constant kinematic viscosity is used.
    if (d_use_constant_kinematic_viscosity_and_ideal_gas_assumptions)
    {
        TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 2);
    }
    // If the constant kinematic dynamic is used.
    else
    {
        TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 1);
    }
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_shear_viscosity = data_shear_viscosity->getGhostBox();
    const hier::IntVector ghostcell_dims_shear_viscosity = ghost_box_shear_viscosity.numberCells();
    
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
    
    hier::IntVector offset_shear_viscosity(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_temperature(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_shear_viscosity = data_shear_viscosity->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_shear_viscosity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_shear_viscosity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_shear_viscosity = num_ghosts_shear_viscosity;
        offset_pressure = num_ghosts_pressure;
        offset_temperature = num_ghosts_temperature;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_shear_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_shear_viscosity = domain.lower() - ghost_box_shear_viscosity.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
    }
    
    // If the constant kinematic viscosity is used.
    if (d_use_constant_kinematic_viscosity_and_ideal_gas_assumptions)
    {
        /*
         * Get the pointers to the cell data.
         */
        
        double* mu = data_shear_viscosity->getPointer(0);
        
        double* p = data_pressure->getPointer(0);
        double* T = data_temperature->getPointer(0);
        
        /*
         * Get the molecular properties.
         */
        
        double nu = *molecular_properties[0];
        double M  = *molecular_properties[1];
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset_0_shear_viscosity = offset_shear_viscosity[0];
            
            const int offset_0_pressure = offset_pressure[0];
            const int offset_0_temperature = offset_temperature[0];
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_shear_viscosity = i + offset_0_shear_viscosity;
                
                const int idx_pressure = i + offset_0_pressure;
                const int idx_temperature = i + offset_0_temperature;
                
                const double R = d_R_u/M;
                
                const double rho = p[idx_pressure]/(R*T[idx_temperature]);
                
                mu[idx_shear_viscosity] = rho*nu;
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
            
            const int offset_0_pressure = offset_pressure[0];
            const int offset_1_pressure = offset_pressure[1];
            const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
            
            const int offset_0_temperature = offset_temperature[0];
            const int offset_1_temperature = offset_temperature[1];
            const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
            
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
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                    
                    const int idx_temperature = (i + offset_0_temperature) +
                        (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                    
                    const double R = d_R_u/M;
                    
                    const double rho = p[idx_pressure]/(R*T[idx_temperature]);
                    
                    mu[idx_shear_viscosity] = rho*nu;
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
                        
                        const int idx_pressure = (i + offset_0_pressure) +
                            (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                            (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                                ghostcell_dim_1_pressure;
                        
                        const int idx_temperature = (i + offset_0_temperature) +
                            (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                            (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                                ghostcell_dim_1_temperature;
                        
                        const double R = d_R_u/M;
                        
                        const double rho = p[idx_pressure]/(R*T[idx_temperature]);
                        
                        mu[idx_shear_viscosity] = rho*nu;
                    }
                }
            }
        }
    }
    // If the constant kinematic dynamic is used.
    else
    {
        /*
         * Get the pointers to the cell data.
         */
        
        double* mu = data_shear_viscosity->getPointer(0);
        
        /*
         * Get the source dynamic viscosity.
         */
        
        double mu_src = *molecular_properties[0];
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and offset.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset_0_shear_viscosity = offset_shear_viscosity[0];
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_shear_viscosity = i + offset_0_shear_viscosity;
                
                mu[idx_shear_viscosity] = mu_src;
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
            
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_shear_viscosity = (i + offset_0_shear_viscosity) +
                        (j + offset_1_shear_viscosity)*ghostcell_dim_0_shear_viscosity;
                    
                    mu[idx_shear_viscosity] = mu_src;
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
                        const int idx_shear_viscosity = (i + offset_0_shear_viscosity) +
                            (j + offset_1_shear_viscosity)*ghostcell_dim_0_shear_viscosity +
                            (k + offset_2_shear_viscosity)*ghostcell_dim_0_shear_viscosity*
                                ghostcell_dim_1_shear_viscosity;
                        
                        mu[idx_shear_viscosity] = mu_src;
                    }
                }
            }
        }
    }
}


/*
 * Compute the shear viscosity.
 */
void
EquationOfShearViscosityConstant::computeShearViscosity(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_shear_viscosity,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_molecular_properties,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_shear_viscosity);
    TBOX_ASSERT(data_molecular_properties);
    TBOX_ASSERT(data_pressure);
    TBOX_ASSERT(data_temperature);
    
    // If the constant kinematic viscosity is used.
    if (d_use_constant_kinematic_viscosity_and_ideal_gas_assumptions)
    {
        TBOX_ASSERT(data_molecular_properties->getDepth() >= 2);
    }
    // If the constant kinematic dynamic is used.
    else
    {
        TBOX_ASSERT(data_molecular_properties->getDepth() >= 1);
    }
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_shear_viscosity = data_shear_viscosity->getGhostBox();
    const hier::IntVector ghostcell_dims_shear_viscosity = ghost_box_shear_viscosity.numberCells();
    
    const hier::Box ghost_box_molecular_properties = data_molecular_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_molecular_properties = ghost_box_molecular_properties.numberCells();
    
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
    
    hier::IntVector offset_shear_viscosity(d_dim);
    hier::IntVector offset_molecular_properties(d_dim);
    hier::IntVector offset_pressure(d_dim);
    hier::IntVector offset_temperature(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_shear_viscosity = data_shear_viscosity->getGhostCellWidth();
        const hier::IntVector num_ghosts_molecular_properties = data_molecular_properties->getGhostCellWidth();
        const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_shear_viscosity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_molecular_properties->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_pressure->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_shear_viscosity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_molecular_properties, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_pressure, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_shear_viscosity = num_ghosts_shear_viscosity;
        offset_molecular_properties = num_ghosts_molecular_properties;
        offset_pressure = num_ghosts_pressure;
        offset_temperature = num_ghosts_temperature;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_shear_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_molecular_properties->getGhostBox().contains(domain));
        TBOX_ASSERT(data_pressure->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_shear_viscosity = domain.lower() - ghost_box_shear_viscosity.lower();
        offset_molecular_properties = domain.lower() - ghost_box_molecular_properties.lower();
        offset_pressure = domain.lower() - ghost_box_pressure.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
    }
    
    // If the constant kinematic viscosity is used.
    if (d_use_constant_kinematic_viscosity_and_ideal_gas_assumptions)
    {
        /*
         * Get the pointers to the cell data.
         */
        
        double* mu = data_shear_viscosity->getPointer(0);
        
        double* nu = data_molecular_properties->getPointer(0);
        double* M  = data_molecular_properties->getPointer(1);
        
        double* p = data_pressure->getPointer(0);
        double* T = data_temperature->getPointer(0);
        
        if (d_dim == tbox::Dimension(1))
        {
            /*
             * Get the local lower index, numbers of cells in each dimension and offsets.
             */
            
            const int domain_lo_0 = domain_lo[0];
            const int domain_dim_0 = domain_dims[0];
            
            const int offset_0_shear_viscosity = offset_shear_viscosity[0];
            const int offset_0_molecular_properties = offset_molecular_properties[0];
            
            const int offset_0_pressure = offset_pressure[0];
            const int offset_0_temperature = offset_temperature[0];
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_shear_viscosity = i + offset_0_shear_viscosity;
                
                const int idx_molecular_properties = i + offset_0_molecular_properties;
                const int idx_pressure = i + offset_0_pressure;
                const int idx_temperature = i + offset_0_temperature;
                
                const double R = d_R_u/M[idx_molecular_properties];
                
                const double rho = p[idx_pressure]/(R*T[idx_temperature]);
                
                mu[idx_shear_viscosity] = rho*nu[idx_molecular_properties];
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
            
            const int offset_0_pressure = offset_pressure[0];
            const int offset_1_pressure = offset_pressure[1];
            const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
            
            const int offset_0_temperature = offset_temperature[0];
            const int offset_1_temperature = offset_temperature[1];
            const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
            
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
                    
                    const int idx_pressure = (i + offset_0_pressure) +
                        (j + offset_1_pressure)*ghostcell_dim_0_pressure;
                    
                    const int idx_temperature = (i + offset_0_temperature) +
                        (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                    
                    const double R = d_R_u/M[idx_molecular_properties];
                    
                    const double rho = p[idx_pressure]/(R*T[idx_temperature]);
                    
                    mu[idx_shear_viscosity] = rho*nu[idx_molecular_properties];
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
                        
                        const int idx_pressure = (i + offset_0_pressure) +
                            (j + offset_1_pressure)*ghostcell_dim_0_pressure +
                            (k + offset_2_pressure)*ghostcell_dim_0_pressure*
                                ghostcell_dim_1_pressure;
                        
                        const int idx_temperature = (i + offset_0_temperature) +
                            (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                            (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                                ghostcell_dim_1_temperature;
                        
                        const double R = d_R_u/M[idx_molecular_properties];
                        
                        const double rho = p[idx_pressure]/(R*T[idx_temperature]);
                        
                        mu[idx_shear_viscosity] = rho*nu[idx_molecular_properties];
                    }
                }
            }
        }
    }
    // If the constant kinematic dynamic is used.
    else
    {
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
}
