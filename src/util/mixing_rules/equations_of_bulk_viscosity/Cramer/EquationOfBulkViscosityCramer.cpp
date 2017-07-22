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
        (gamma - 1.0)*(gamma - 1.0)*c_v_v*A_v*exp(B_v/(pow(T, 1.0/3.0)) + C_v/(pow(T, 2.0/3.0)));
    
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
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_bulk_viscosity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_bulk_viscosity = data_bulk_viscosity->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_bulk_viscosity =
        data_bulk_viscosity->getGhostBox().numberCells();
    
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
        
        num_ghosts_min = num_ghosts_bulk_viscosity;
        num_ghosts_min = hier::IntVector::min(num_ghosts_temperature, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_bulk_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
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
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_bulk_viscosity = num_ghosts_bulk_viscosity[0];
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_bulk_viscosity = i + num_ghosts_0_bulk_viscosity;
            const int idx_temperature = i + num_ghosts_0_temperature;
            
            mu_v[idx_bulk_viscosity] = A_r + B_r*T[idx_temperature] + 
                (gamma - 1.0)*(gamma - 1.0)*c_v_v*A_v*exp(B_v/(pow(T[idx_temperature], 1.0/3.0)) +
                    C_v/(pow(T[idx_temperature], 2.0/3.0)));
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
        
        const int num_ghosts_0_bulk_viscosity = num_ghosts_bulk_viscosity[0];
        const int num_ghosts_1_bulk_viscosity = num_ghosts_bulk_viscosity[1];
        const int ghostcell_dim_0_bulk_viscosity = ghostcell_dims_bulk_viscosity[0];
        
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
                const int idx_bulk_viscosity = (i + num_ghosts_0_bulk_viscosity) +
                    (j + num_ghosts_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity;
                
                const int idx_temperature = (i + num_ghosts_0_temperature) +
                    (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                
                mu_v[idx_bulk_viscosity] = A_r + B_r*T[idx_temperature] + 
                    (gamma - 1.0)*(gamma - 1.0)*c_v_v*A_v*exp(B_v/(pow(T[idx_temperature], 1.0/3.0)) +
                        C_v/(pow(T[idx_temperature], 2.0/3.0)));
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
        
        const int num_ghosts_0_bulk_viscosity = num_ghosts_bulk_viscosity[0];
        const int num_ghosts_1_bulk_viscosity = num_ghosts_bulk_viscosity[1];
        const int num_ghosts_2_bulk_viscosity = num_ghosts_bulk_viscosity[2];
        const int ghostcell_dim_0_bulk_viscosity = ghostcell_dims_bulk_viscosity[0];
        const int ghostcell_dim_1_bulk_viscosity = ghostcell_dims_bulk_viscosity[1];
        
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
                    const int idx_bulk_viscosity = (i + num_ghosts_0_bulk_viscosity) +
                        (j + num_ghosts_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity +
                        (k + num_ghosts_2_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity*
                            ghostcell_dim_1_bulk_viscosity;
                    
                    const int idx_temperature = (i + num_ghosts_0_temperature) +
                        (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                        (k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    mu_v[idx_bulk_viscosity] = A_r + B_r*T[idx_temperature] + 
                        (gamma - 1.0)*(gamma - 1.0)*c_v_v*A_v*exp(B_v/(pow(T[idx_temperature], 1.0/3.0)) +
                            C_v/(pow(T[idx_temperature], 2.0/3.0)));
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
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_bulk_viscosity->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_temperature->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(data_molecular_properties->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_bulk_viscosity = data_bulk_viscosity->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_bulk_viscosity =
        data_bulk_viscosity->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_temperature =
        data_temperature->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_molecular_properties = data_molecular_properties->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_molecular_properties =
        data_molecular_properties->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_bulk_viscosity;
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
        TBOX_ASSERT(data_bulk_viscosity->getGhostBox().contains(domain));
        TBOX_ASSERT(data_temperature->getGhostBox().contains(domain));
        TBOX_ASSERT(data_molecular_properties->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
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
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_bulk_viscosity = num_ghosts_bulk_viscosity[0];
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_0_molecular_properties = num_ghosts_molecular_properties[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
        {
            // Compute the linear indices.
            const int idx_bulk_viscosity = i + num_ghosts_0_bulk_viscosity;
            const int idx_temperature = i + num_ghosts_0_temperature;
            const int idx_molecular_properties = i + num_ghosts_0_molecular_properties;
            
            mu_v[idx_bulk_viscosity] = A_r[idx_molecular_properties] +
                B_r[idx_molecular_properties]*T[idx_temperature] + 
                (gamma[idx_molecular_properties] - 1.0)*(gamma[idx_molecular_properties] - 1.0)*
                    c_v_v[idx_molecular_properties]*A_v[idx_molecular_properties]*
                        exp(B_v[idx_molecular_properties]/(pow(T[idx_temperature], 1.0/3.0)) +
                            C_v[idx_molecular_properties]/(pow(T[idx_temperature], 2.0/3.0)));
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
        
        const int num_ghosts_0_bulk_viscosity = num_ghosts_bulk_viscosity[0];
        const int num_ghosts_1_bulk_viscosity = num_ghosts_bulk_viscosity[1];
        const int ghostcell_dim_0_bulk_viscosity = ghostcell_dims_bulk_viscosity[0];
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        
        const int num_ghosts_0_molecular_properties = num_ghosts_molecular_properties[0];
        const int num_ghosts_1_molecular_properties = num_ghosts_molecular_properties[1];
        const int ghostcell_dim_0_molecular_properties = ghostcell_dims_molecular_properties[0];
        
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_bulk_viscosity = (i + num_ghosts_0_bulk_viscosity) +
                    (j + num_ghosts_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity;
                
                const int idx_temperature = (i + num_ghosts_0_temperature) +
                    (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                
                const int idx_molecular_properties = (i + num_ghosts_0_molecular_properties) +
                    (j + num_ghosts_1_molecular_properties)*ghostcell_dim_0_molecular_properties;
                
                mu_v[idx_bulk_viscosity] = A_r[idx_molecular_properties] +
                    B_r[idx_molecular_properties]*T[idx_temperature] + 
                    (gamma[idx_molecular_properties] - 1.0)*(gamma[idx_molecular_properties] - 1.0)*
                        c_v_v[idx_molecular_properties]*A_v[idx_molecular_properties]*
                            exp(B_v[idx_molecular_properties]/(pow(T[idx_temperature], 1.0/3.0)) +
                                C_v[idx_molecular_properties]/(pow(T[idx_temperature], 2.0/3.0)));
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
        
        const int num_ghosts_0_bulk_viscosity = num_ghosts_bulk_viscosity[0];
        const int num_ghosts_1_bulk_viscosity = num_ghosts_bulk_viscosity[1];
        const int num_ghosts_2_bulk_viscosity = num_ghosts_bulk_viscosity[2];
        const int ghostcell_dim_0_bulk_viscosity = ghostcell_dims_bulk_viscosity[0];
        const int ghostcell_dim_1_bulk_viscosity = ghostcell_dims_bulk_viscosity[1];
        
        const int num_ghosts_0_temperature = num_ghosts_temperature[0];
        const int num_ghosts_1_temperature = num_ghosts_temperature[1];
        const int num_ghosts_2_temperature = num_ghosts_temperature[2];
        const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
        const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
        
        const int num_ghosts_0_molecular_properties = num_ghosts_molecular_properties[0];
        const int num_ghosts_1_molecular_properties = num_ghosts_molecular_properties[1];
        const int num_ghosts_2_molecular_properties = num_ghosts_molecular_properties[2];
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
                    const int idx_bulk_viscosity = (i + num_ghosts_0_bulk_viscosity) +
                        (j + num_ghosts_1_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity +
                        (k + num_ghosts_2_bulk_viscosity)*ghostcell_dim_0_bulk_viscosity*
                            ghostcell_dim_1_bulk_viscosity;
                    
                    const int idx_temperature = (i + num_ghosts_0_temperature) +
                        (j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                        (k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                            ghostcell_dim_1_temperature;
                    
                    const int idx_molecular_properties = (i + num_ghosts_0_molecular_properties) +
                        (j + num_ghosts_1_molecular_properties)*ghostcell_dim_0_molecular_properties +
                        (k + num_ghosts_2_molecular_properties)*ghostcell_dim_0_molecular_properties*
                            ghostcell_dim_1_molecular_properties;
                    
                    mu_v[idx_bulk_viscosity] = A_r[idx_molecular_properties] +
                        B_r[idx_molecular_properties]*T[idx_temperature] + 
                        (gamma[idx_molecular_properties] - 1.0)*(gamma[idx_molecular_properties] - 1.0)*
                            c_v_v[idx_molecular_properties]*A_v[idx_molecular_properties]*
                                exp(B_v[idx_molecular_properties]/(pow(T[idx_temperature], 1.0/3.0)) +
                                    C_v[idx_molecular_properties]/(pow(T[idx_temperature], 2.0/3.0)));
                }
            }
        }
    }
}
