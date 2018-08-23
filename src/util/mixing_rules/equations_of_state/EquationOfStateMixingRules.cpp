#include "util/mixing_rules/equations_of_state/EquationOfStateMixingRules.hpp"

/*
 * Get the molecular weight of a species.
 */
double
EquationOfStateMixingRules::getSpeciesMolecularWeight(
    const int species_index) const
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(species_index >= 0);
    TBOX_ASSERT(species_index < d_num_species);
#endif
    
    return d_species_M[species_index];
}


/*
 * Compute the molecular weight of mixture.
 */
double
EquationOfStateMixingRules::getMixtureMolecularWeight(
    const std::vector<const double*>& mass_fractions) const
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    if (static_cast<int>(mass_fractions.size()) != d_num_species)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Number of mass fractions provided is not"
            << "equal to the total number of species."
            << std::endl);
    }
#endif
    
    double M = double(0);
    
    for (int si = 0; si < d_num_species; si++)
    {
        const double& Y = *(mass_fractions[si]);
        
        M += Y/d_species_M[si];
    }
    
    return double(1)/M; 
}


/*
 * Compute the molecular weight of mixture.
 */
void
EquationOfStateMixingRules::computeMixtureMolecularWeight(
    boost::shared_ptr<pdat::CellData<double> >& data_mixture_molecular_weight,
    const boost::shared_ptr<pdat::CellData<double> >& data_mass_fractions,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_mixture_molecular_weight);
    TBOX_ASSERT(data_mass_fractions);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_mixture_molecular_weight->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_mixture_molecular_weight = data_mixture_molecular_weight->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_mixture_molecular_weight =
        data_mixture_molecular_weight->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_mass_fractions =
        data_mass_fractions->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_molecular_weight;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_mixture_molecular_weight->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* M = data_mixture_molecular_weight->getPointer(0);
    
    std::vector<const double*> Y;
    Y.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Y.push_back(data_mass_fractions->getPointer(si));
    }
    
    /*
     * Fill zeros for rho.
     */
    
    if (domain.empty())
    {
        data_mixture_molecular_weight->fillAll(double(0));
    }
    else
    {
        data_mixture_molecular_weight->fillAll(double(0), domain);
    }
    
    computeMixtureMolecularWeight(
        M,
        Y,
        num_ghosts_mixture_molecular_weight,
        num_ghosts_mass_fractions,
        ghostcell_dims_mixture_molecular_weight,
        ghostcell_dims_mass_fractions,
        domain_lo,
        domain_dims);
}


/*
 * Compute the molecular weight of mixture.
 */
void
EquationOfStateMixingRules::computeMixtureMolecularWeight(
    boost::shared_ptr<pdat::SideData<double> >& data_mixture_molecular_weight,
    const boost::shared_ptr<pdat::SideData<double> >& data_mass_fractions,
    int side_normal,
    const hier::Box& domain) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_mixture_molecular_weight);
    TBOX_ASSERT(data_mass_fractions);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = data_mixture_molecular_weight->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_mass_fractions->getBox().numberCells() == interior_dims);
#endif
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_mixture_molecular_weight = data_mixture_molecular_weight->getGhostCellWidth();
    hier::IntVector ghostcell_dims_mixture_molecular_weight =
        data_mixture_molecular_weight->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
    hier::IntVector ghostcell_dims_mass_fractions =
        data_mass_fractions->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = num_ghosts_mixture_molecular_weight;
        num_ghosts_min = hier::IntVector::min(num_ghosts_mass_fractions, num_ghosts_min);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_mixture_molecular_weight->getGhostBox().contains(domain));
        TBOX_ASSERT(data_mass_fractions->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_mixture_molecular_weight->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_mass_fractions->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_mixture_molecular_weight[side_normal]++;
    ghostcell_dims_mass_fractions[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* M = data_mixture_molecular_weight->getPointer(side_normal, 0);
    
    std::vector<const double*> Y;
    Y.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Y.push_back(data_mass_fractions->getPointer(side_normal, si));
    }
    
    /*
     * Fill zeros for rho.
     */
    
    if (domain.empty())
    {
        data_mixture_molecular_weight->fillAll(double(0));
    }
    else
    {
        data_mixture_molecular_weight->fillAll(double(0), domain);
    }
    
    computeMixtureMolecularWeight(
        M,
        Y,
        num_ghosts_mixture_molecular_weight,
        num_ghosts_mass_fractions,
        ghostcell_dims_mixture_molecular_weight,
        ghostcell_dims_mass_fractions,
        domain_lo,
        domain_dims);
}


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
    
    std::vector<const double*> Z_rho;
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
    
    computeMixtureDensity(
        rho,
        Z_rho,
        num_ghosts_mixture_density,
        num_ghosts_partial_densities,
        ghostcell_dims_mixture_density,
        ghostcell_dims_partial_densities,
        domain_lo,
        domain_dims);
}


/*
 * Helper function to compute the density of mixture given the partial densities.
 */
void
EquationOfStateMixingRules::computeMixtureDensity(
    boost::shared_ptr<pdat::SideData<double> >& data_mixture_density,
    const boost::shared_ptr<pdat::SideData<double> >& data_partial_densities,
    int side_normal,
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
    hier::IntVector ghostcell_dims_mixture_density =
        data_mixture_density->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_partial_densities = data_partial_densities->getGhostCellWidth();
    hier::IntVector ghostcell_dims_partial_densities =
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
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(side_normal < d_dim.getValue());
    
    TBOX_ASSERT(data_mixture_density->getDirectionVector()[side_normal] > 0);
    TBOX_ASSERT(data_partial_densities->getDirectionVector()[side_normal] > 0);
#endif
    
    ghostcell_dims_mixture_density[side_normal]++;
    ghostcell_dims_partial_densities[side_normal]++;
    domain_dims[side_normal]++;
    
    /*
     * Get the pointers to the cell data.
     */
    
    double* rho = data_mixture_density->getPointer(side_normal, 0);
    
    std::vector<const double*> Z_rho;
    Z_rho.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Z_rho.push_back(data_partial_densities->getPointer(side_normal, si));
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
    
    computeMixtureDensity(
        rho,
        Z_rho,
        num_ghosts_mixture_density,
        num_ghosts_partial_densities,
        ghostcell_dims_mixture_density,
        ghostcell_dims_partial_densities,
        domain_lo,
        domain_dims);
}


/*
 * Compute the molecular weight of mixture given the mass fractions.
 */
void
EquationOfStateMixingRules::computeMixtureMolecularWeight(
    double* const M,
    const std::vector<const double*>& Y,
    const hier::IntVector& num_ghosts_mixture_molecular_weight,
    const hier::IntVector& num_ghosts_mass_fractions,
    const hier::IntVector& ghostcell_dims_mixture_molecular_weight,
    const hier::IntVector& ghostcell_dims_mass_fractions,
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
        
        const int num_ghosts_0_mixture_molecular_weight = num_ghosts_mixture_molecular_weight[0];
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        
        for (int si = 0; si < d_num_species; si++)
        {
            const double M_i = d_species_M[si];
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0;
                 i < domain_lo_0 + domain_dim_0;
                 i++)
            {
                // Compute the linear indices.
                const int idx_mixture_molecular_weight = i + num_ghosts_0_mixture_molecular_weight;
                const int idx_mass_fractions = i + num_ghosts_0_mass_fractions;
                
                M[idx_mixture_molecular_weight] += Y[si][idx_mass_fractions]/M_i;
            }
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0;
             i < domain_lo_0 + domain_dim_0;
             i++)
        {
            // Compute the linear index.
            const int idx_mixture_molecular_weight = i + num_ghosts_0_mixture_molecular_weight;
            
            M[idx_mixture_molecular_weight] = double(1)/M[idx_mixture_molecular_weight];
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
        
        const int num_ghosts_0_mixture_molecular_weight = num_ghosts_mixture_molecular_weight[0];
        const int num_ghosts_1_mixture_molecular_weight = num_ghosts_mixture_molecular_weight[1];
        const int ghostcell_dim_0_mixture_molecular_weight = ghostcell_dims_mixture_molecular_weight[0];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        
        for (int si = 0; si < d_num_species; si++)
        {
            const double M_i = d_species_M[si];
            
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
                    const int idx_mixture_molecular_weight = (i + num_ghosts_0_mixture_molecular_weight) +
                        (j + num_ghosts_1_mixture_molecular_weight)*ghostcell_dim_0_mixture_molecular_weight;
                    
                    const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                        (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                    
                    M[idx_mixture_molecular_weight] += Y[si][idx_mass_fractions]/M_i;
                }
            }
        }
        
            
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
                // Compute the linear index.
                const int idx_mixture_molecular_weight = (i + num_ghosts_0_mixture_molecular_weight) +
                    (j + num_ghosts_1_mixture_molecular_weight)*ghostcell_dim_0_mixture_molecular_weight;
                
                M[idx_mixture_molecular_weight] = double(1)/M[idx_mixture_molecular_weight];
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
        
        const int num_ghosts_0_mixture_molecular_weight = num_ghosts_mixture_molecular_weight[0];
        const int num_ghosts_1_mixture_molecular_weight = num_ghosts_mixture_molecular_weight[1];
        const int num_ghosts_2_mixture_molecular_weight = num_ghosts_mixture_molecular_weight[2];
        const int ghostcell_dim_0_mixture_molecular_weight = ghostcell_dims_mixture_molecular_weight[0];
        const int ghostcell_dim_1_mixture_molecular_weight = ghostcell_dims_mixture_molecular_weight[1];
        
        const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
        const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
        const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
        const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
        const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
        
        for (int si = 0; si < d_num_species; si++)
        {
            const double M_i = d_species_M[si];
            
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
                        const int idx_mixture_molecular_weight = (i + num_ghosts_0_mixture_molecular_weight) +
                            (j + num_ghosts_1_mixture_molecular_weight)*ghostcell_dim_0_mixture_molecular_weight +
                            (k + num_ghosts_2_mixture_molecular_weight)*ghostcell_dim_0_mixture_molecular_weight*
                                ghostcell_dim_1_mixture_molecular_weight;
                        
                        const int idx_mass_fractions = (i + num_ghosts_0_mass_fractions) +
                            (j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                            (k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                ghostcell_dim_1_mass_fractions;
                        
                        M[idx_mixture_molecular_weight] += Y[si][idx_mass_fractions]/M_i;
                    }
                }
            }
        }
        
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
                    // Compute the linear index.
                    const int idx_mixture_molecular_weight = (i + num_ghosts_0_mixture_molecular_weight) +
                        (j + num_ghosts_1_mixture_molecular_weight)*ghostcell_dim_0_mixture_molecular_weight +
                        (k + num_ghosts_2_mixture_molecular_weight)*ghostcell_dim_0_mixture_molecular_weight*
                            ghostcell_dim_1_mixture_molecular_weight;
                    
                    M[idx_mixture_molecular_weight] = double(1)/M[idx_mixture_molecular_weight];
                }
            }
        }
    }
}


/*
 * Compute the density of mixture given the partial densities.
 */
void
EquationOfStateMixingRules::computeMixtureDensity(
    double* const rho,
    const std::vector<const double*>& Z_rho,
    const hier::IntVector& num_ghosts_mixture_density,
    const hier::IntVector& num_ghosts_partial_densities,
    const hier::IntVector& ghostcell_dims_mixture_density,
    const hier::IntVector& ghostcell_dims_partial_densities,
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

