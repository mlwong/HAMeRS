#include "util/mixing_rules/equations_of_state/EquationOfStateMixingRules.hpp"

/*
 * Helper function to compute the density of mixture given the partial densities.
 */
double
EquationOfStateMixingRules::getMixtureDensity(
    const std::vector<const double*>& partial_densities) const
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    if (static_cast<int>(partial_densities.size()) != d_num_species)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of partial densities provided is not"
            << "equal to the total number of species."
            << std::endl);
    }
#endif
    
    double rho = double(0);
    
    for (int si = 0; si < d_num_species; si++)
    {
        const double& Z_rho = *(partial_densities[si]);
        
        rho += Z_rho;
    }
    
    return rho;            
}


/*
 * Helper function to compute the density of mixture given the partial densities.
 */
void
EquationOfStateMixingRules::computeMixtureDensity(
    boost::shared_ptr<pdat::CellData<double> >& data_mixture_density,
    const boost::shared_ptr<pdat::CellData<double> >& data_partial_densities,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_mixture_density);
    TBOX_ASSERT(data_partial_densities);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_mixture_density->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_partial_densities->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_mixture_density = data_mixture_density->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_mixture_density =
        data_mixture_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_partial_densities = data_partial_densities->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_partial_densities =
        data_partial_densities->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_density;
        num_ghosts_min = hier::IntVector::min(num_ghosts_partial_densities, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_mixture_density->getGhostBox().contains(domain));
        TBOX_ASSERT(data_partial_densities->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* rho = data_mixture_density->getPointer(0);
    
    std::vector<double*> Z_rho;
    Z_rho.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Z_rho.push_back(data_partial_densities->getPointer(si));
    }
    
    /*
     * Fill zeros for rho.
     */
    
    if (domain.empty())
    {
        data_mixture_density->fillAll(double(0));
    }
    else
    {
        data_mixture_density->fillAll(double(0), domain);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_mixture_density = num_ghosts_mixture_density[0];
        const int num_ghosts_0_partial_densities = num_ghosts_partial_densities[0];
        
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0;
                 i < domain_lo_0 + domain_dim_0;
                 i++)
            {
                // Compute the linear indices.
                const int idx_mixture_density = i + num_ghosts_0_mixture_density;
                const int idx_partial_densities = i + num_ghosts_0_partial_densities;
                
                rho[idx_mixture_density] += Z_rho[si][idx_partial_densities];
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
        
        const int num_ghosts_0_mixture_density = num_ghosts_mixture_density[0];
        const int num_ghosts_1_mixture_density = num_ghosts_mixture_density[1];
        const int ghostcell_dim_0_mixture_density = ghostcell_dims_mixture_density[0];
        
        const int num_ghosts_0_partial_densities = num_ghosts_partial_densities[0];
        const int num_ghosts_1_partial_densities = num_ghosts_partial_densities[1];
        const int ghostcell_dim_0_partial_densities = ghostcell_dims_partial_densities[0];
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = domain_lo_1;
                 j < domain_lo_1 + domain_dim_1;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = domain_lo_0;
                     i < domain_lo_0 + domain_dim_0;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_mixture_density = (i + num_ghosts_0_mixture_density) +
                        (j + num_ghosts_1_mixture_density)*ghostcell_dim_0_mixture_density;
                    
                    const int idx_partial_densities = (i + num_ghosts_0_partial_densities) +
                        (j + num_ghosts_1_partial_densities)*ghostcell_dim_0_partial_densities;
                    
                    rho[idx_mixture_density] += Z_rho[si][idx_partial_densities];
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
        
        const int num_ghosts_0_mixture_density = num_ghosts_mixture_density[0];
        const int num_ghosts_1_mixture_density = num_ghosts_mixture_density[1];
        const int num_ghosts_2_mixture_density = num_ghosts_mixture_density[2];
        const int ghostcell_dim_0_mixture_density = ghostcell_dims_mixture_density[0];
        const int ghostcell_dim_1_mixture_density = ghostcell_dims_mixture_density[1];
        
        const int num_ghosts_0_partial_densities = num_ghosts_partial_densities[0];
        const int num_ghosts_1_partial_densities = num_ghosts_partial_densities[1];
        const int num_ghosts_2_partial_densities = num_ghosts_partial_densities[2];
        const int ghostcell_dim_0_partial_densities = ghostcell_dims_partial_densities[0];
        const int ghostcell_dim_1_partial_densities = ghostcell_dims_partial_densities[1];
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = domain_lo_2;
                 k < domain_lo_2 + domain_dim_2;
                 k++)
            {
                for (int j = domain_lo_1;
                     j < domain_lo_1 + domain_dim_1;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = domain_lo_0;
                         i < domain_lo_0 + domain_dim_0;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_mixture_density = (i + num_ghosts_0_mixture_density) +
                            (j + num_ghosts_1_mixture_density)*ghostcell_dim_0_mixture_density +
                            (k + num_ghosts_2_mixture_density)*ghostcell_dim_0_mixture_density*
                                ghostcell_dim_1_mixture_density;
                        
                        const int idx_partial_densities = (i + num_ghosts_0_partial_densities) +
                            (j + num_ghosts_1_partial_densities)*ghostcell_dim_0_partial_densities +
                            (k + num_ghosts_2_partial_densities)*ghostcell_dim_0_partial_densities*
                                ghostcell_dim_1_partial_densities;
                        
                        rho[idx_mixture_density] += Z_rho[si][idx_partial_densities];
                    }
                }
            }
        }
    }
}

