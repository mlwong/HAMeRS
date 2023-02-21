#include "util/mixing_rules/equations_of_shear_viscosity/Chapman-Enskog/EquationOfShearViscosityChapmanEnskog.hpp"

#include <cmath>

/*
 * Print all characteristics of the equation of shear viscosity class.
 */
void
EquationOfShearViscosityChapmanEnskog::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfShearViscosityChapmanEnskog object..."
       << std::endl;
    
    os << std::endl;
    
    os << "EquationOfShearViscosityChapmanEnskog: this = "
       << (EquationOfShearViscosityChapmanEnskog *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the shear viscosity.
 */
double
EquationOfShearViscosityChapmanEnskog::getShearViscosity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& molecular_properties) const
{
    NULL_USE(pressure);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 3);
#endif
    
    double mu = double(0);
    
    const double& epsilon_by_k = *(molecular_properties[0]);
    const double& sigma = *(molecular_properties[1]);
    const double& M = *(molecular_properties[2]);
    
    const double& T = *temperature;
    
    const double A = double(1.16145);
    const double B = double(-0.14874);
    const double C = double(0.52487);
    const double D = double(-0.7732);
    const double E = double(2.16178);
    const double F = double(-2.43787);
    
    const double T_star = T/epsilon_by_k;
    const double Omega = A*pow(T_star, B) + C*exp(D*T_star) + E*exp(F*T_star);
    
    mu = double(2.6693e-6)*sqrt(M*T)/(Omega*sigma*sigma);
    
    return mu;
}


/*
 * Compute the shear viscosity.
 */
void
EquationOfShearViscosityChapmanEnskog::computeShearViscosity(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_shear_viscosity,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
    const std::vector<const double*>& molecular_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_shear_viscosity);
    TBOX_ASSERT(data_temperature);
    
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 3);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_shear_viscosity = data_shear_viscosity->getGhostBox();
    const hier::IntVector ghostcell_dims_shear_viscosity = ghost_box_shear_viscosity.numberCells();
    
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    const hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_shear_viscosity(d_dim);
    hier::IntVector offset_temperature(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_shear_viscosity = data_shear_viscosity->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_shear_viscosity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_shear_viscosity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_shear_viscosity = num_ghosts_shear_viscosity;
        offset_temperature = num_ghosts_temperature;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_shear_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_shear_viscosity = domain.lower() - ghost_box_shear_viscosity.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* mu = data_shear_viscosity->getPointer(0);
    double* T = data_temperature->getPointer(0);
    
    const double& epsilon_by_k = *(molecular_properties[0]);
    const double& sigma = *(molecular_properties[1]);
    const double& M = *(molecular_properties[2]);
    
    const double A = double(1.16145);
    const double B = double(-0.14874);
    const double C = double(0.52487);
    const double D = double(-0.7732);
    const double E = double(2.16178);
    const double F = double(-2.43787);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_shear_viscosity = offset_shear_viscosity[0];
        const int offset_0_temperature = offset_temperature[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_shear_viscosity = i + offset_0_shear_viscosity;
            const int idx_temperature = i + offset_0_temperature;
            
            const double T_star = T[idx_temperature]/epsilon_by_k;
            const double Omega = A*pow(T_star, B) + C*exp(D*T_star) + E*exp(F*T_star);
            
            mu[idx_shear_viscosity] = double(2.6693e-6)*sqrt(M*T[idx_temperature])/(Omega*sigma*sigma);
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
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_shear_viscosity = (i + offset_0_shear_viscosity) +
                    (j + offset_1_shear_viscosity)*ghostcell_dim_0_shear_viscosity;
                
                const int idx_temperature = (i + offset_0_temperature) +
                    (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                
                const double T_star = T[idx_temperature]/epsilon_by_k;
                const double Omega = A*pow(T_star, B) + C*exp(D*T_star) + E*exp(F*T_star);
                
                mu[idx_shear_viscosity] = double(2.6693e-6)*sqrt(M*T[idx_temperature])/(Omega*sigma*sigma);
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
                    const int idx_shear_viscosity = (i + offset_0_shear_viscosity) +
                        (j + offset_1_shear_viscosity)*ghostcell_dim_0_shear_viscosity +
                        (k + offset_2_shear_viscosity)*ghostcell_dim_0_shear_viscosity*
                            ghostcell_dim_1_shear_viscosity;
                    
                    const int idx_temperature = (i + offset_0_temperature) +
                        (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                        (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    const double T_star = T[idx_temperature]/epsilon_by_k;
                    const double Omega = A*pow(T_star, B) + C*exp(D*T_star) + E*exp(F*T_star);
                    
                    mu[idx_shear_viscosity] = double(2.6693e-6)*sqrt(M*T[idx_temperature])/(Omega*sigma*sigma);
                }
            }
        }
    }
}


/*
 * Compute the shear viscosity.
 */
void
EquationOfShearViscosityChapmanEnskog::computeShearViscosity(
    HAMERS_SHARED_PTR<pdat::CellData<double> >& data_shear_viscosity,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_pressure,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_temperature,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_molecular_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_shear_viscosity);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_molecular_properties);
    
    TBOX_ASSERT(data_molecular_properties->getDepth() >= 3);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_shear_viscosity = data_shear_viscosity->getGhostBox();
    const hier::IntVector ghostcell_dims_shear_viscosity = ghost_box_shear_viscosity.numberCells();
    
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    const hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    const hier::Box ghost_box_molecular_properties = data_molecular_properties->getGhostBox();
    const hier::IntVector ghostcell_dims_molecular_properties = ghost_box_molecular_properties.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_shear_viscosity(d_dim);
    hier::IntVector offset_temperature(d_dim);
    hier::IntVector offset_molecular_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_shear_viscosity = data_shear_viscosity->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_molecular_properties = data_molecular_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_shear_viscosity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_molecular_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_shear_viscosity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_molecular_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_shear_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_molecular_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_shear_viscosity = domain.lower() - ghost_box_shear_viscosity.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
        offset_molecular_properties = domain.lower() - ghost_box_molecular_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* mu = data_shear_viscosity->getPointer(0);
    double* T = data_temperature->getPointer(0);
    
    double* epsilon_by_k = data_molecular_properties->getPointer(0);
    double* sigma = data_molecular_properties->getPointer(1);
    double* M = data_molecular_properties->getPointer(2);
    
    const double A = double(1.16145);
    const double B = double(-0.14874);
    const double C = double(0.52487);
    const double D = double(-0.7732);
    const double E = double(2.16178);
    const double F = double(-2.43787);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_shear_viscosity = offset_shear_viscosity[0];
        const int offset_0_temperature = offset_temperature[0];
        const int offset_0_molecular_properties = offset_molecular_properties[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_shear_viscosity = i + offset_0_shear_viscosity;
            const int idx_temperature = i + offset_0_temperature;
            const int idx_molecular_properties = i + offset_0_molecular_properties;
            
            const double T_star = T[idx_temperature]/epsilon_by_k[idx_molecular_properties];
            const double Omega = A*pow(T_star, B) + C*exp(D*T_star) + E*exp(F*T_star);
            
            mu[idx_shear_viscosity] = double(2.6693e-6)*sqrt(M[idx_molecular_properties]*T[idx_temperature])/
                (Omega*sigma[idx_molecular_properties]*sigma[idx_molecular_properties]);
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
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        const int offset_0_molecular_properties = offset_molecular_properties[0];
        const int offset_1_molecular_properties = offset_molecular_properties[1];
        const int ghostcell_dim_0_molecular_properties = ghostcell_dims_molecular_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_shear_viscosity = (i + offset_0_shear_viscosity) +
                    (j + offset_1_shear_viscosity)*ghostcell_dim_0_shear_viscosity;
                
                const int idx_temperature = (i + offset_0_temperature) +
                    (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                
                const int idx_molecular_properties = (i + offset_0_molecular_properties) +
                    (j + offset_1_molecular_properties)*ghostcell_dim_0_molecular_properties;
                
                const double T_star = T[idx_temperature]/epsilon_by_k[idx_molecular_properties];
                const double Omega = A*pow(T_star, B) + C*exp(D*T_star) + E*exp(F*T_star);
                
                mu[idx_shear_viscosity] = double(2.6693e-6)*sqrt(M[idx_molecular_properties]*T[idx_temperature])/
                    (Omega*sigma[idx_molecular_properties]*sigma[idx_molecular_properties]);
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
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int offset_2_temperature = offset_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
        const int offset_0_molecular_properties = offset_molecular_properties[0];
        const int offset_1_molecular_properties = offset_molecular_properties[1];
        const int offset_2_molecular_properties = offset_molecular_properties[2];
        const int ghostcell_dim_0_molecular_properties = ghostcell_dims_molecular_properties[0];
        const int ghostcell_dim_1_molecular_properties = ghostcell_dims_molecular_properties[1];
        
        for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
        {
            for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_shear_viscosity = (i + offset_0_shear_viscosity) +
                        (j + offset_1_shear_viscosity)*ghostcell_dim_0_shear_viscosity +
                        (k + offset_2_shear_viscosity)*ghostcell_dim_0_shear_viscosity*
                            ghostcell_dim_1_shear_viscosity;
                    
                    const int idx_temperature = (i + offset_0_temperature) +
                        (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                        (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    const int idx_molecular_properties = (i + offset_0_molecular_properties) +
                        (j + offset_1_molecular_properties)*ghostcell_dim_0_molecular_properties +
                        (k + offset_2_molecular_properties)*ghostcell_dim_0_molecular_properties*
                            ghostcell_dim_1_molecular_properties;
                    
                    const double T_star = T[idx_temperature]/epsilon_by_k[idx_molecular_properties];
                    const double Omega = A*pow(T_star, B) + C*exp(D*T_star) + E*exp(F*T_star);
                    
                    mu[idx_shear_viscosity] = double(2.6693e-6)*sqrt(M[idx_molecular_properties]*T[idx_temperature])/
                        (Omega*sigma[idx_molecular_properties]*sigma[idx_molecular_properties]);
                }
            }
        }
    }
}
