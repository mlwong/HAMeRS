#include "util/mixing_rules/equations_of_bulk_viscosity/Cramer/EquationOfBulkViscosityCramer.hpp"

#include <cmath>

/*
 * Print all characteristics of the equation of bulk viscosity class.
 */
void
EquationOfBulkViscosityCramer::printClassData(
    std::ostream& os) const
{
    os << "\nPrint EquationOfBulkViscosityCramer object..."
       << std::endl;
       
    os << std::endl;
    
    os << "EquationOfBulkViscosityCramer: this = "
       << (EquationOfBulkViscosityCramer *)this
       << std::endl;
    
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Compute the bulk viscosity.
 */
double
EquationOfBulkViscosityCramer::getBulkViscosity(
    const double* const pressure,
    const double* const temperature,
    const std::vector<const double*>& molecular_properties) const
{
    NULL_USE(pressure);
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 7);
#endif
    
    const double& T = *temperature;
    
    const double& gamma = *(molecular_properties[0]);
    
    const double& A_r = *(molecular_properties[1]);
    const double& B_r = *(molecular_properties[2]);
    
    const double& c_v_v = *(molecular_properties[3]);
    const double& A_v   = *(molecular_properties[4]);
    const double& B_v   = *(molecular_properties[5]);
    const double& C_v   = *(molecular_properties[6]);
    
    double mu_v = A_r + B_r*T +
        (gamma - double(1))*(gamma - double(1))*c_v_v*A_v*exp(B_v/(pow(T, double(1)/double(3))) +
        C_v/(pow(T, double(2)/double(3))));
    
    return mu_v;
}


/*
 * Compute the bulk viscosity.
 */
void
EquationOfBulkViscosityCramer::computeBulkViscosity(
    boost::shared_ptr<pdat::CellData<double> >& data_bulk_viscosity,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const std::vector<const double*>& molecular_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_bulk_viscosity);
    TBOX_ASSERT(data_temperature);
    
    TBOX_ASSERT(static_cast<int>(molecular_properties.size()) >= 7);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_bulk_viscosity = data_bulk_viscosity->getGhostBox();
    const hier::IntVector ghostcell_dims_bulk_viscosity = ghost_box_bulk_viscosity.numberCells();
    
    const hier::Box ghost_box_temperature = data_temperature->getGhostBox();
    const hier::IntVector ghostcell_dims_temperature = ghost_box_temperature.numberCells();
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     * Also, get the offsets.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    hier::IntVector offset_bulk_viscosity(d_dim);
    hier::IntVector offset_temperature(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_bulk_viscosity = data_bulk_viscosity->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_bulk_viscosity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_bulk_viscosity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_bulk_viscosity = num_ghosts_bulk_viscosity;
        offset_temperature = num_ghosts_temperature;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_bulk_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_bulk_viscosity = domain.lower() - ghost_box_bulk_viscosity.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* mu_v = data_bulk_viscosity->getPointer(0);
    double* T = data_temperature->getPointer(0);
    
    const double& gamma = *(molecular_properties[0]);
    
    const double& A_r = *(molecular_properties[1]);
    const double& B_r = *(molecular_properties[2]);
    
    const double& c_v_v = *(molecular_properties[3]);
    const double& A_v   = *(molecular_properties[4]);
    const double& B_v   = *(molecular_properties[5]);
    const double& C_v   = *(molecular_properties[6]);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_bulk_viscosity = offset_bulk_viscosity[0];
        const int offset_0_temperature = offset_temperature[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_bulk_viscosity = i + offset_0_bulk_viscosity;
            const int idx_temperature = i + offset_0_temperature;
            
            mu_v[idx_bulk_viscosity] = A_r + B_r*T[idx_temperature] + 
                (gamma - double(1))*(gamma - double(1))*c_v_v*A_v*
                exp(B_v/(pow(T[idx_temperature], double(1)/double(3))) +
                C_v/(pow(T[idx_temperature], double(2)/double(3))));
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
                const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                    (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity;
                
                const int idx_temperature = (i + offset_0_temperature) +
                    (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                
                mu_v[idx_bulk_viscosity] = A_r + B_r*T[idx_temperature] + 
                    (gamma - double(1))*(gamma - double(1))*c_v_v*A_v*
                    exp(B_v/(pow(T[idx_temperature], double(1)/double(3))) +
                    C_v/(pow(T[idx_temperature], double(2)/double(3))));
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
                    const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                        (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity +
                        (k + offset_2_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity*
                            ghostcell_dim_1_bulk_viscosity;
                    
                    const int idx_temperature = (i + offset_0_temperature) +
                        (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                        (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    mu_v[idx_bulk_viscosity] = A_r + B_r*T[idx_temperature] + 
                        (gamma - double(1))*(gamma - double(1))*c_v_v*
                        A_v*exp(B_v/(pow(T[idx_temperature], double(1)/double(3))) +
                        C_v/(pow(T[idx_temperature], double(2)/double(3))));
                }
            }
        }
    }
}


/*
 * Compute the bulk viscosity.
 */
void
EquationOfBulkViscosityCramer::computeBulkViscosity(
    boost::shared_ptr<pdat::CellData<double> >& data_bulk_viscosity,
    const boost::shared_ptr<pdat::CellData<double> >& data_pressure,
    const boost::shared_ptr<pdat::CellData<double> >& data_temperature,
    const boost::shared_ptr<pdat::CellData<double> >& data_molecular_properties,
    const hier::Box& domain) const
{
    NULL_USE(data_pressure);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_bulk_viscosity);
    TBOX_ASSERT(data_temperature);
    TBOX_ASSERT(data_molecular_properties);
    
    TBOX_ASSERT(data_molecular_properties->getDepth() >= 7);
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_bulk_viscosity = data_bulk_viscosity->getGhostBox();
    const hier::IntVector ghostcell_dims_bulk_viscosity = ghost_box_bulk_viscosity.numberCells();
    
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
    
    hier::IntVector offset_bulk_viscosity(d_dim);
    hier::IntVector offset_temperature(d_dim);
    hier::IntVector offset_molecular_properties(d_dim);
    
    if (domain.empty())
    {
        // Get the numbers of ghost cells.
        const hier::IntVector num_ghosts_bulk_viscosity = data_bulk_viscosity->getGhostCellWidth();
        const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
        const hier::IntVector num_ghosts_molecular_properties = data_molecular_properties->getGhostCellWidth();
        
        // Get the box that covers the interior of patch.
        const hier::Box interior_box = data_bulk_viscosity->getBox();
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_temperature->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_molecular_properties->getBox().isSpatiallyEqual(interior_box));
#endif
        
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_bulk_viscosity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        num_ghosts_min = hier::IntVector::min(num_ghosts_molecular_properties, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
        
        offset_bulk_viscosity = num_ghosts_bulk_viscosity;
        offset_temperature = num_ghosts_temperature;
        offset_molecular_properties = num_ghosts_molecular_properties;
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_bulk_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_molecular_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = hier::IntVector::getZero(d_dim);
        domain_dims = domain.numberCells();
        
        offset_bulk_viscosity = domain.lower() - ghost_box_bulk_viscosity.lower();
        offset_temperature = domain.lower() - ghost_box_temperature.lower();
        offset_molecular_properties = domain.lower() - ghost_box_molecular_properties.lower();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* mu_v = data_bulk_viscosity->getPointer(0);
    double* T = data_temperature->getPointer(0);
    
    double* gamma = data_molecular_properties->getPointer(0);
    
    double* A_r = data_molecular_properties->getPointer(1);
    double* B_r = data_molecular_properties->getPointer(2);
    
    double* c_v_v = data_molecular_properties->getPointer(3);
    double* A_v   = data_molecular_properties->getPointer(4);
    double* B_v   = data_molecular_properties->getPointer(5);
    double* C_v   = data_molecular_properties->getPointer(6);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and offsets.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int offset_0_bulk_viscosity = offset_bulk_viscosity[0];
        const int offset_0_temperature = offset_temperature[0];
        const int offset_0_molecular_properties = offset_molecular_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_bulk_viscosity = i + offset_0_bulk_viscosity;
            const int idx_temperature = i + offset_0_temperature;
            const int idx_molecular_properties = i + offset_0_molecular_properties;
            
            mu_v[idx_bulk_viscosity] = A_r[idx_molecular_properties] +
                B_r[idx_molecular_properties]*T[idx_temperature] + 
                (gamma[idx_molecular_properties] - double(1))*(gamma[idx_molecular_properties] - double(1))*
                c_v_v[idx_molecular_properties]*A_v[idx_molecular_properties]*
                exp(B_v[idx_molecular_properties]/(pow(T[idx_temperature], double(1)/double(3))) +
                C_v[idx_molecular_properties]/(pow(T[idx_temperature], double(2)/double(3))));
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
        
        const int offset_0_temperature = offset_temperature[0];
        const int offset_1_temperature = offset_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
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
                const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                    (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity;
                
                const int idx_temperature = (i + offset_0_temperature) +
                    (j + offset_1_temperature)*ghostcell_dim_0_temperature;
                
                const int idx_molecular_properties = (i + offset_0_molecular_properties) +
                    (j + offset_1_molecular_properties)*ghostcell_dim_0_molecular_properties;
                
                mu_v[idx_bulk_viscosity] = A_r[idx_molecular_properties] +
                    B_r[idx_molecular_properties]*T[idx_temperature] + 
                    (gamma[idx_molecular_properties] - double(1))*(gamma[idx_molecular_properties] - double(1))*
                    c_v_v[idx_molecular_properties]*A_v[idx_molecular_properties]*
                    exp(B_v[idx_molecular_properties]/(pow(T[idx_temperature], double(1)/double(3))) +
                    C_v[idx_molecular_properties]/(pow(T[idx_temperature], double(2)/double(3))));
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
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_bulk_viscosity = (i + offset_0_bulk_viscosity) +
                        (j + offset_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity +
                        (k + offset_2_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity*
                            ghostcell_dim_1_bulk_viscosity;
                    
                    const int idx_temperature = (i + offset_0_temperature) +
                        (j + offset_1_temperature)*ghostcell_dim_0_temperature +
                        (k + offset_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    const int idx_molecular_properties = (i + offset_0_molecular_properties) +
                        (j + offset_1_molecular_properties)*ghostcell_dim_0_molecular_properties +
                        (k + offset_2_molecular_properties)*ghostcell_dim_0_molecular_properties*
                            ghostcell_dim_1_molecular_properties;
                    
                    mu_v[idx_bulk_viscosity] = A_r[idx_molecular_properties] +
                        B_r[idx_molecular_properties]*T[idx_temperature] + 
                        (gamma[idx_molecular_properties] - double(1))*(gamma[idx_molecular_properties] - double(1))*
                        c_v_v[idx_molecular_properties]*A_v[idx_molecular_properties]*
                        exp(B_v[idx_molecular_properties]/(pow(T[idx_temperature], double(1)/double(3))) +
                        C_v[idx_molecular_properties]/(pow(T[idx_temperature], double(2)/double(3))));
                }
            }
        }
    }
}
