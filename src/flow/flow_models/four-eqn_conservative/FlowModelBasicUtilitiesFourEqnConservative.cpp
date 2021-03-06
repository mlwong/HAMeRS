#include "flow/flow_models/four-eqn_conservative/FlowModelBasicUtilitiesFourEqnConservative.hpp"

/*
 * Convert conservative variables to primitive variables.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::convertConservativeVariablesToPrimitiveVariables(
    const std::vector<const double*>& conservative_variables,
    const std::vector<double*>& primitive_variables)
{
    const std::vector<const double*>& Q = conservative_variables;
    const std::vector<double*>&       V = primitive_variables;
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    if (static_cast<int>(Q.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "Number of elements in conservative variables is not correct."
            << std::endl);
    }
    
    if (static_cast<int>(V.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "Number of elements in primitive variables is not correct."
            << std::endl);
    }
#endif
    
    // Compute the mixture density.
    std::vector<const double*> rho_Y_ptr;
    rho_Y_ptr.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        rho_Y_ptr.push_back(Q[si]);
    }
    const double rho = d_equation_of_state_mixing_rules->getMixtureDensity(
        rho_Y_ptr);
    
    // Compute the specific internal energy.
    double epsilon = double(0);
    if (d_dim == tbox::Dimension(1))
    {
        epsilon = ((*Q[d_num_species + d_dim.getValue()]) -
            double(1)/double(2)*((*Q[d_num_species])*(*Q[d_num_species]))/rho)/rho;
    }
    else if (d_dim == tbox::Dimension(2))
    {
        epsilon = ((*Q[d_num_species + d_dim.getValue()]) -
            double(1)/double(2)*((*Q[d_num_species])*(*Q[d_num_species]) +
            (*Q[d_num_species + 1])*(*Q[d_num_species + 1]))/rho)/rho;
    }
    else if (d_dim == tbox::Dimension(3))
    {
        epsilon = ((*Q[d_num_species + d_dim.getValue()]) -
            double(1)/double(2)*((*Q[d_num_species])*(*Q[d_num_species]) +
            (*Q[d_num_species + 1])*(*Q[d_num_species + 1]) +
            (*Q[d_num_species + 2])*(*Q[d_num_species + 2]))/rho)/rho;
    }
    
    // Compute the mass fractions.
    double Y[d_num_species];
    for (int si = 0; si < d_num_species; si++)
    {
        Y[si] = (*Q[si])/rho;
    }
    
    // Get the pointers to the mass fractions.
    std::vector<const double*> Y_ptr;
    Y_ptr.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Y_ptr.push_back(&Y[si]);
    }
    
    // Compute the pressure.
    const double p = d_equation_of_state_mixing_rules->getPressure(
        &rho,
        &epsilon,
        Y_ptr);
    
    // Convert the conservative variables to primitive variables.
    for (int si = 0; si < d_num_species; si++)
    {
        *V[si] = *Q[si];
    }
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        *V[d_num_species + di] = (*Q[d_num_species + di])/rho;
    }
    *V[d_num_species + d_dim.getValue()] = p;
}


/*
 * Convert conservative variables to primitive variables.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::convertConservativeVariablesToPrimitiveVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    
    /*
     * Get the dimensions of box that covers the interior of patch.
     */
    
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_primitive_var = primitive_variables[0]->
        getGhostCellWidth();
    
    const hier::IntVector num_ghosts_conservative_var = conservative_variables[0]->
        getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of of the variables.
     */
    
    const hier::IntVector ghostcell_dims_primitive_var = primitive_variables[0]->
        getGhostBox().numberCells();
    
    const hier::IntVector ghostcell_dims_conservative_var = conservative_variables[0]->
        getGhostBox().numberCells();
    
    /*
     * Get the size of variables.
     */
    
    int d_num_eqn_primitive_var = 0;
    int d_num_eqn_conservative_var = 0;
    
    /*
     * Check the size of variables.
     */
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        d_num_eqn_primitive_var += primitive_variables[vi]->getDepth();
    }
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        d_num_eqn_conservative_var += conservative_variables[vi]->getDepth();
    }
    
    if (d_num_eqn_primitive_var != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    
    if (d_num_eqn_conservative_var != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_primitive_var =
            primitive_variables[vi]->getBox().numberCells();
        
        if (interior_dims_primitive_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "convertConservativeVariablesToPrimitiveVariables()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_conservative_var =
            conservative_variables[vi]->getBox().numberCells();
        
        if (interior_dims_conservative_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "convertConservativeVariablesToPrimitiveVariables()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int vi = 1; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        if (num_ghosts_primitive_var != primitive_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "convertConservativeVariablesToPrimitiveVariables()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    for (int vi = 1; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        if (num_ghosts_conservative_var != conservative_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "convertConservativeVariablesToPrimitiveVariables()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_primitive_var > num_ghosts_conservative_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "The ghost cell width of primitive variables is larger than that of conservative variables."
            << std::endl);
    }
    
    /*
     * Declare the pointers to the primitive variables and conservative variables.
     */
    
    std::vector<double*> V;
    V.resize(d_num_eqn_primitive_var);
    
    std::vector<double*> Q;
    Q.resize(d_num_eqn_conservative_var);
    
    int count_eqn = 0;
    
    /*
     * Convert conservative variables to primitive variables.
     */
    
    // Create the temporary side data.
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_density(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_internal_energy(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_pressure(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mass_fractions(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_var));
    
    data_density->fillAll(double(0));
    
    double* rho     = nullptr;
    double* epsilon = nullptr;
    double* p       = nullptr;
    
    std::vector<double*> Y;
    Y.resize(d_num_species);
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_primitive_var    = num_ghosts_primitive_var[0];
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        
        /*
         * Convert conservative variables to primitive variables in the x-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(0, 0);
        epsilon = data_internal_energy->getPointer(0, 0);
        p       = data_pressure->getPointer(0, 0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(0, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_conservative_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_conservative_var = i + num_ghosts_0_conservative_var;
                
                rho[idx_conservative_var] += Q[si][idx_conservative_var];
            }
        }
        
        // Compute the internal energy.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear index.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            epsilon[idx_conservative_var] = (Q[d_num_species + d_dim.getValue()][idx_conservative_var] -
                double(1)/double(2)*(Q[d_num_species][idx_conservative_var]*Q[d_num_species][idx_conservative_var])/
                rho[idx_conservative_var])/rho[idx_conservative_var];
        }
        
        // Compute the mass fractions.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_conservative_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                    i++)
            {
                // Compute the linear index.
                const int idx_conservative_var = i + num_ghosts_0_conservative_var;
                
                Y[si][idx_conservative_var] = Q[si][idx_conservative_var]/rho[idx_conservative_var];
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            data_mass_fractions,
            0);
        
        // Set the partial densities.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                 i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
                const int idx_conservative_var = i + num_ghosts_0_conservative_var;
                
                V[si][idx_primitive_var] = Q[si][idx_conservative_var];
            }
        }
        
        // Set the velocity.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            V[d_num_species][idx_primitive_var] = Q[d_num_species][idx_conservative_var]/rho[idx_conservative_var];
        }
        
        // Set the pressure.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            V[d_num_species + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        const int num_ghosts_1_primitive_var = num_ghosts_primitive_var[1];
        const int ghostcell_dim_0_primitive_var = ghostcell_dims_primitive_var[0];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[1];
        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[0];
        
        /*
         * Convert conservative variables to primitive variables in the x-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(0, 0);
        epsilon = data_internal_energy->getPointer(0, 0);
        p       = data_pressure->getPointer(0, 0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(0, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_conservative_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                    
                    rho[idx_conservative_var] += Q[si][idx_conservative_var];
                }
            }
        }
        
        // Compute the internal energy.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_conservative_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                    i++)
            {
                // Compute the linear index.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                epsilon[idx_conservative_var] = (Q[d_num_species + d_dim.getValue()][idx_conservative_var] -
                    double(1)/double(2)*(Q[d_num_species][idx_conservative_var]*Q[d_num_species][idx_conservative_var] +
                    Q[d_num_species + 1][idx_conservative_var]*Q[d_num_species + 1][idx_conservative_var])/
                    rho[idx_conservative_var])/rho[idx_conservative_var];
            }
        }
        
        // Compute the mass fractions.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_conservative_var;
                        i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                        i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                    
                    Y[si][idx_conservative_var] = Q[si][idx_conservative_var]/rho[idx_conservative_var];
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            data_mass_fractions,
            0);
        
        // Set the partial densities.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_primitive_var;
                        i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                        i++)
                {
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                    
                    V[si][idx_primitive_var] = Q[si][idx_conservative_var];
                }
            }
        }
        
        // Set the velocity.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                V[d_num_species][idx_primitive_var] = Q[d_num_species][idx_conservative_var]/rho[idx_conservative_var];
                V[d_num_species + 1][idx_primitive_var] = Q[d_num_species + 1][idx_conservative_var]/rho[idx_conservative_var];
            }
        }
        
        // Set the pressure.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                V[d_num_species + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
            }
        }
        
        /*
         * Convert conservative variables to primitive variables in the y-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(1, 0);
        epsilon = data_internal_energy->getPointer(1, 0);
        p       = data_pressure->getPointer(1, 0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(1, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                    
                    rho[idx_conservative_var] += Q[si][idx_conservative_var];
                }
            }
        }
        
        // Compute the internal energy.
        for (int j = -num_ghosts_1_conservative_var;
                j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                epsilon[idx_conservative_var] = (Q[d_num_species + d_dim.getValue()][idx_conservative_var] -
                    double(1)/double(2)*(Q[d_num_species][idx_conservative_var]*Q[d_num_species][idx_conservative_var] +
                    Q[d_num_species + 1][idx_conservative_var]*Q[d_num_species + 1][idx_conservative_var])/
                    rho[idx_conservative_var])/rho[idx_conservative_var];
            }
        }
        
        // Compute the mass fractions.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                    
                    Y[si][idx_conservative_var] = Q[si][idx_conservative_var]/rho[idx_conservative_var];
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            data_mass_fractions,
            1);
        
        // Set the partial densities.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                    
                    V[si][idx_primitive_var] = Q[si][idx_conservative_var];
                }
            }
        }
        
        // Set the velocity.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                V[d_num_species][idx_primitive_var] = Q[d_num_species][idx_conservative_var]/rho[idx_conservative_var];
                V[d_num_species + 1][idx_primitive_var] = Q[d_num_species + 1][idx_conservative_var]/rho[idx_conservative_var];
            }
        }
        
        // Set the pressure.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                V[d_num_species + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        const int num_ghosts_1_primitive_var = num_ghosts_primitive_var[1];
        const int num_ghosts_2_primitive_var = num_ghosts_primitive_var[2];
        const int ghostcell_dim_0_primitive_var = ghostcell_dims_primitive_var[0];
        const int ghostcell_dim_1_primitive_var = ghostcell_dims_primitive_var[1];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[1];
        const int num_ghosts_2_conservative_var = num_ghosts_conservative_var[2];
        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[0];
        const int ghostcell_dim_1_conservative_var = ghostcell_dims_conservative_var[1];
        
        /*
         * Convert conservative variables to primitive variables in the x-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(0, 0);
        epsilon = data_internal_energy->getPointer(0, 0);
        p       = data_pressure->getPointer(0, 0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(0, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_conservative_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                            (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                                ghostcell_dim_1_conservative_var;
                        
                        rho[idx_conservative_var] += Q[si][idx_conservative_var];
                    }
                }
            }
        }
        
        // Compute the internal energy.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_conservative_var;
                        i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                        i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    epsilon[idx_conservative_var] = (Q[d_num_species + d_dim.getValue()][idx_conservative_var] -
                        double(1)/double(2)*(Q[d_num_species][idx_conservative_var]*Q[d_num_species][idx_conservative_var] +
                        Q[d_num_species + 1][idx_conservative_var]*Q[d_num_species + 1][idx_conservative_var] +
                        Q[d_num_species + 2][idx_conservative_var]*Q[d_num_species + 2][idx_conservative_var])/
                        rho[idx_conservative_var])/rho[idx_conservative_var];
                }
            }
        }
        
        // Compute the mass fractions.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_conservative_var;
                            i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                            i++)
                    {
                        // Compute the linear index.
                        const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                            (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                                ghostcell_dim_1_conservative_var;
                        
                        Y[si][idx_conservative_var] = Q[si][idx_conservative_var]/rho[idx_conservative_var];
                    }
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            data_mass_fractions,
            0);
        
        // Set the partial densities.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_primitive_var;
                            i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                            i++)
                    {
                        // Compute the linear indices.
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                            (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                                ghostcell_dim_1_primitive_var;
                        
                        const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                            (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                                ghostcell_dim_1_conservative_var;
                        
                        V[si][idx_primitive_var] = Q[si][idx_conservative_var];
                    }
                }
            }
        }
        
        // Set the velocity.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_primitive_var;
                        i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                        i++)
                {
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    V[d_num_species][idx_primitive_var] = Q[d_num_species][idx_conservative_var]/rho[idx_conservative_var];
                    V[d_num_species + 1][idx_primitive_var] = Q[d_num_species + 1][idx_conservative_var]/rho[idx_conservative_var];
                    V[d_num_species + 2][idx_primitive_var] = Q[d_num_species + 2][idx_conservative_var]/rho[idx_conservative_var];
                }
            }
        }
        
        // Set the pressure.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_primitive_var;
                        i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                        i++)
                {
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    V[d_num_species + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
                }
            }
        }
        
        /*
         * Convert conservative variables to primitive variables in the y-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(1, 0);
        epsilon = data_internal_energy->getPointer(1, 0);
        p       = data_pressure->getPointer(1, 0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(1, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_conservative_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                (ghostcell_dim_1_conservative_var + 1);
                        
                        rho[idx_conservative_var] += Q[si][idx_conservative_var];
                    }
                }
            }
        }
        
        // Compute the internal energy.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                    j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                    j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    epsilon[idx_conservative_var] = (Q[d_num_species + d_dim.getValue()][idx_conservative_var] -
                        double(1)/double(2)*(Q[d_num_species][idx_conservative_var]*Q[d_num_species][idx_conservative_var] +
                        Q[d_num_species + 1][idx_conservative_var]*Q[d_num_species + 1][idx_conservative_var] +
                        Q[d_num_species + 2][idx_conservative_var]*Q[d_num_species + 2][idx_conservative_var])/
                        rho[idx_conservative_var])/rho[idx_conservative_var];
                }
            }
        }
        
        // Compute the mass fractions.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_conservative_var;
                        j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                        j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                (ghostcell_dim_1_conservative_var + 1);
                        
                        Y[si][idx_conservative_var] = Q[si][idx_conservative_var]/rho[idx_conservative_var];
                    }
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            data_mass_fractions,
            1);
        
        // Set the partial densities.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_primitive_var;
                        j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                        j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                (ghostcell_dim_1_primitive_var + 1);
                        
                        const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                (ghostcell_dim_1_conservative_var + 1);
                        
                        V[si][idx_primitive_var] = Q[si][idx_conservative_var];
                    }
                }
            }
        }
        
        // Set the velocity.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                    j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                    j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    V[d_num_species][idx_primitive_var] = Q[d_num_species][idx_conservative_var]/rho[idx_conservative_var];
                    V[d_num_species + 1][idx_primitive_var] = Q[d_num_species + 1][idx_conservative_var]/rho[idx_conservative_var];
                    V[d_num_species + 2][idx_primitive_var] = Q[d_num_species + 2][idx_conservative_var]/rho[idx_conservative_var];
                }
            }
        }
        
        // Set the pressure.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                    j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                    j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    V[d_num_species + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
                }
            }
        }
        
        /*
         * Convert conservative variables to primitive variables in the z-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(2, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(2, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(2, 0);
        epsilon = data_internal_energy->getPointer(2, 0);
        p       = data_pressure->getPointer(2, 0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(2, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_conservative_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                ghostcell_dim_1_conservative_var;
                        
                        rho[idx_conservative_var] += Q[si][idx_conservative_var];
                    }
                }
            }
        }
        
        // Compute the internal energy.
        for (int k = -num_ghosts_2_conservative_var;
                k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
                k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    epsilon[idx_conservative_var] = (Q[d_num_species + d_dim.getValue()][idx_conservative_var] -
                        double(1)/double(2)*(Q[d_num_species][idx_conservative_var]*Q[d_num_species][idx_conservative_var] +
                        Q[d_num_species + 1][idx_conservative_var]*Q[d_num_species + 1][idx_conservative_var] +
                        Q[d_num_species + 2][idx_conservative_var]*Q[d_num_species + 2][idx_conservative_var])/
                        rho[idx_conservative_var])/rho[idx_conservative_var];
                }
            }
        }
        
        // Compute the mass fractions.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_conservative_var;
                    k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
                    k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                ghostcell_dim_1_conservative_var;
                        
                        Y[si][idx_conservative_var] = Q[si][idx_conservative_var]/rho[idx_conservative_var];
                    }
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            data_mass_fractions,
            2);
        
        // Set the partial densities.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_primitive_var;
                    k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
                    k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                ghostcell_dim_1_primitive_var;
                        
                        const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                ghostcell_dim_1_conservative_var;
                        
                        V[si][idx_primitive_var] = Q[si][idx_conservative_var];
                    }
                }
            }
        }
        
        // Set the velocity.
        for (int k = -num_ghosts_2_primitive_var;
                k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
                k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    V[d_num_species][idx_primitive_var] = Q[d_num_species][idx_conservative_var]/rho[idx_conservative_var];
                    V[d_num_species + 1][idx_primitive_var] = Q[d_num_species + 1][idx_conservative_var]/rho[idx_conservative_var];
                    V[d_num_species + 2][idx_primitive_var] = Q[d_num_species + 2][idx_conservative_var]/rho[idx_conservative_var];
                }
            }
        }
        
        // Set the pressure.
        for (int k = -num_ghosts_2_primitive_var;
                k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
                k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    V[d_num_species + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
                }
            }
        }
    }
}


/*
 * Convert primitive variables to conservative variables.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::convertPrimitiveVariablesToConservativeVariables(
    const std::vector<const double*>& primitive_variables,
    const std::vector<double*>& conservative_variables)
{
    const std::vector<const double*>& V = primitive_variables;
    const std::vector<double*>&       Q = conservative_variables;
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    if (static_cast<int>(V.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "Number of elements in primitive variables is not correct."
            << std::endl);
    }
    
    if (static_cast<int>(Q.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "Number of elements in conservative variables is not correct."
            << std::endl);
    }
#endif
    
    // Compute the mixture density.
    std::vector<const double*> rho_Y_ptr;
    rho_Y_ptr.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        rho_Y_ptr.push_back(V[si]);
    }
    const double rho = d_equation_of_state_mixing_rules->getMixtureDensity(
        rho_Y_ptr);
    
    // Compute the mass fractions.
    double Y[d_num_species];
    for (int si = 0; si < d_num_species; si++)
    {
        Y[si] = (*V[si])/rho;
    }
    
    // Get the pointers to the mass fractions.
    std::vector<const double*> Y_ptr;
    Y_ptr.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        Y_ptr.push_back(&Y[si]);
    }
    
    // Compute the total energy.
    const double epsilon = d_equation_of_state_mixing_rules->getInternalEnergy(
        &rho,
        V[d_num_species + d_dim.getValue()],
        Y_ptr);
    
    double E = double(0);
    if (d_dim == tbox::Dimension(1))
    {
        E = rho*(epsilon + double(1)/double(2)*((*V[d_num_species])*(*V[d_num_species])));
    }
    else if (d_dim == tbox::Dimension(2))
    {
        E = rho*(epsilon + double(1)/double(2)*(
            (*V[d_num_species])*(*V[d_num_species]) +
            (*V[d_num_species + 1])*(*V[d_num_species + 1])));
    }
    else if (d_dim == tbox::Dimension(3))
    {
        E = rho*(epsilon + double(1)/double(2)*(
            (*V[d_num_species])*(*V[d_num_species]) +
            (*V[d_num_species + 1])*(*V[d_num_species + 1]) + 
            (*V[d_num_species + 2])*(*V[d_num_species + 2])));
    }
    
    // Convert the primitive variables to conservative variables.
    for (int si = 0; si < d_num_species; si++)
    {
        *Q[si] = *V[si];
    }
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        *Q[d_num_species + di] = rho*(*V[d_num_species + di]);
    }
    *Q[d_num_species + d_dim.getValue()] = E;
}


/*
 * Convert primitive variables to conservative variables.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::convertPrimitiveVariablesToConservativeVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    
    /*
     * Get the dimensions of box that covers the interior of patch.
     */
    
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_conservative_var = conservative_variables[0]->
        getGhostCellWidth();
    
    const hier::IntVector num_ghosts_primitive_var = primitive_variables[0]->
        getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of of the variables.
     */
    
    const hier::IntVector ghostcell_dims_conservative_var = conservative_variables[0]->
        getGhostBox().numberCells();
    
    const hier::IntVector ghostcell_dims_primitive_var = primitive_variables[0]->
        getGhostBox().numberCells();
    
    /*
     * Get the size of variables.
     */
    
    int d_num_eqn_conservative_var = 0;
    int d_num_eqn_primitive_var = 0;
    
    /*
     * Check the size of variables.
     */
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        d_num_eqn_conservative_var += conservative_variables[vi]->getDepth();
    }
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        d_num_eqn_primitive_var += primitive_variables[vi]->getDepth();
    }
    
    if (d_num_eqn_conservative_var != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    
    if (d_num_eqn_primitive_var != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_conservative_var =
            conservative_variables[vi]->getBox().numberCells();
        
        if (interior_dims_conservative_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "convertPrimitiveVariablesToConservativeVariables()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_primitive_var =
            primitive_variables[vi]->getBox().numberCells();
        
        if (interior_dims_primitive_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "convertPrimitiveVariablesToConservativeVariables()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int vi = 1; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        if (num_ghosts_conservative_var != conservative_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "convertPrimitiveVariablesToConservativeVariables()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    for (int vi = 1; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        if (num_ghosts_primitive_var != primitive_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "convertPrimitiveVariablesToConservativeVariables()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_conservative_var > num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "The ghost cell width of conservative variables is larger than that of primitive variables."
            << std::endl);
    }
    
    /*
     * Declare the pointers to the conservative variables and primitive variables.
     */
    
    std::vector<double*> Q;
    Q.resize(d_num_eqn_conservative_var);
    
    std::vector<double*> V;
    V.resize(d_num_eqn_primitive_var);
    
    int count_eqn = 0;
    
    /*
     * Convert primitive variables to conservative variables.
     */
    
    // Create the temporary side data.
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_density(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_pressure(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_internal_energy(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_mass_fractions(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_var));
    
    data_density->fillAll(double(0));
    
    double* rho     = nullptr;
    double* p       = nullptr;
    double* epsilon = nullptr;
    
    std::vector<double*> Y;
    Y.resize(d_num_species);
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_0_primitive_var    = num_ghosts_primitive_var[0];
        
        /*
         * Convert primitive variables to conservative variables in the x-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(0, 0);
        epsilon = data_internal_energy->getPointer(0, 0);
        p       = data_pressure->getPointer(0, 0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(0, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_primitive_var = i + num_ghosts_0_primitive_var;
                
                rho[idx_primitive_var] += V[si][idx_primitive_var];
            }
        }
        
        // Get the pressure.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear index.
            const int idx_primitive_var = i + num_ghosts_0_primitive_var;
            
            p[idx_primitive_var] = V[d_num_species + d_dim.getValue()][idx_primitive_var];
        }
        
        // Compute the mass fractions.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                    i++)
            {
                // Compute the linear index.
                const int idx_primitive_var = i + num_ghosts_0_primitive_var;
                
                Y[si][idx_primitive_var] = V[si][idx_primitive_var]/rho[idx_primitive_var];
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            data_mass_fractions,
            0);
        
        // Set the partial densities.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_conservative_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                 i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = i + num_ghosts_0_conservative_var;
                const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
                
                Q[si][idx_conservative_var] = V[si][idx_primitive_var];
            }
        }
        
        // Set the momentum.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            
            Q[d_num_species][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species][idx_primitive_var];
        }
        
        // Set the total energy.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            
            Q[d_num_species + d_dim.getValue()][idx_conservative_var] = rho[idx_primitive_var]*
                (epsilon[idx_primitive_var] + double(1)/double(2)*
                V[d_num_species][idx_primitive_var]*V[d_num_species][idx_primitive_var]);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[1];
        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[0];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        const int num_ghosts_1_primitive_var = num_ghosts_primitive_var[1];
        const int ghostcell_dim_0_primitive_var = ghostcell_dims_primitive_var[0];
        
        /*
         * Convert primitive variables to conservative variables in the x-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(0, 0);
        epsilon = data_internal_energy->getPointer(0, 0);
        p       = data_pressure->getPointer(0, 0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(0, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_primitive_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                    
                    rho[idx_primitive_var] += V[si][idx_primitive_var];
                }
            }
        }
        
        // Get the pressure.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                    i++)
            {
                // Compute the linear index.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                p[idx_primitive_var] = V[d_num_species + d_dim.getValue()][idx_primitive_var];
            }
        }
        
        // Compute the mass fractions.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_primitive_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                    
                    Y[si][idx_primitive_var] = V[si][idx_primitive_var]/rho[idx_primitive_var];
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            data_mass_fractions,
            0);
        
        // Set the partial densities.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_conservative_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                    
                    Q[si][idx_conservative_var] = V[si][idx_primitive_var];
                }
            }
        }
        
        // Set the momentum.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_conservative_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                Q[d_num_species][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species][idx_primitive_var];
                Q[d_num_species + 1][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species + 1][idx_primitive_var];
            }
        }
        
        // Set the total energy.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_conservative_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                Q[d_num_species + d_dim.getValue()][idx_conservative_var] = rho[idx_primitive_var]*
                    (epsilon[idx_primitive_var] + double(1)/double(2)*(
                    V[d_num_species][idx_primitive_var]*V[d_num_species][idx_primitive_var] +
                    V[d_num_species + 1][idx_primitive_var]*V[d_num_species + 1][idx_primitive_var]));
            }
        }
        
        /*
         * Convert primitive variables to conservative variables in the y-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(1, 0);
        epsilon = data_internal_energy->getPointer(1, 0);
        p       = data_pressure->getPointer(1, 0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(1, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                    
                    rho[idx_primitive_var] += V[si][idx_primitive_var];
                }
            }
        }
        
        // Get the pressure.
        for (int j = -num_ghosts_1_primitive_var;
                j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                p[idx_primitive_var] = V[d_num_species + d_dim.getValue()][idx_primitive_var];
            }
        }
        
        // Compute the mass fractions.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                    
                    Y[si][idx_primitive_var] = V[si][idx_primitive_var]/rho[idx_primitive_var];
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            data_mass_fractions,
            1);
        
        // Set the partial densities.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                    
                    Q[si][idx_conservative_var] = V[si][idx_primitive_var];
                }
            }
        }
        
        // Set the momentum.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                Q[d_num_species][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species][idx_primitive_var];
                Q[d_num_species + 1][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species + 1][idx_primitive_var];
            }
        }
        
        // Set the total energy.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                Q[d_num_species + d_dim.getValue()][idx_conservative_var] = rho[idx_primitive_var]*
                    (epsilon[idx_primitive_var] + double(1)/double(2)*(
                    V[d_num_species][idx_primitive_var]*V[d_num_species][idx_primitive_var] +
                    V[d_num_species + 1][idx_primitive_var]*V[d_num_species + 1][idx_primitive_var]));
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[1];
        const int num_ghosts_2_conservative_var = num_ghosts_conservative_var[2];
        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[0];
        const int ghostcell_dim_1_conservative_var = ghostcell_dims_conservative_var[1];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        const int num_ghosts_1_primitive_var = num_ghosts_primitive_var[1];
        const int num_ghosts_2_primitive_var = num_ghosts_primitive_var[2];
        const int ghostcell_dim_0_primitive_var = ghostcell_dims_primitive_var[0];
        const int ghostcell_dim_1_primitive_var = ghostcell_dims_primitive_var[1];
        
        /*
         * Convert primitive variables to conservative variables in the x-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(0, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(0, 0);
        epsilon = data_internal_energy->getPointer(0, 0);
        p       = data_pressure->getPointer(0, 0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(0, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_primitive_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                            (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                                ghostcell_dim_1_primitive_var;
                        
                        rho[idx_primitive_var] += V[si][idx_primitive_var];
                    }
                }
            }
        }
        
        // Get the pressure.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_primitive_var;
                        i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                        i++)
                {
                    // Compute the linear index.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    p[idx_primitive_var] = V[d_num_species + d_dim.getValue()][idx_primitive_var];
                }
            }
        }
        
        // Compute the mass fractions.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_primitive_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                            (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                                ghostcell_dim_1_primitive_var;
                        
                        Y[si][idx_primitive_var] = V[si][idx_primitive_var]/rho[idx_primitive_var];
                    }
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            data_mass_fractions,
            0);
        
        // Set the partial densities.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_conservative_var;
                            i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                            i++)
                    {
                        // Compute the linear indices.
                        const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                            (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                                ghostcell_dim_1_conservative_var;
                        
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                            (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                                ghostcell_dim_1_primitive_var;
                        
                        Q[si][idx_conservative_var] = V[si][idx_primitive_var];
                    }
                }
            }
        }
        
        // Set the momentum.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_conservative_var;
                        i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                        i++)
                {
                    // Compute the linear indices.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    Q[d_num_species][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species][idx_primitive_var];
                    Q[d_num_species + 1][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species + 1][idx_primitive_var];
                    Q[d_num_species + 2][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species + 2][idx_primitive_var];
                }
            }
        }
        
        // Set the total energy.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_conservative_var;
                        i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                        i++)
                {
                    // Compute the linear indices.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    Q[d_num_species + d_dim.getValue()][idx_conservative_var] = rho[idx_primitive_var]*
                        (epsilon[idx_primitive_var] + double(1)/double(2)*(
                        V[d_num_species][idx_primitive_var]*V[d_num_species][idx_primitive_var] + 
                        V[d_num_species + 1][idx_primitive_var]*V[d_num_species + 1][idx_primitive_var] +
                        V[d_num_species + 2][idx_primitive_var]*V[d_num_species + 2][idx_primitive_var]));
                }
            }
        }
        
        /*
         * Convert primitive variables to conservative variables in the y-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(1, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(1, 0);
        epsilon = data_internal_energy->getPointer(1, 0);
        p       = data_pressure->getPointer(1, 0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(1, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_primitive_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                (ghostcell_dim_1_primitive_var + 1);
                        
                        rho[idx_primitive_var] += V[si][idx_primitive_var];
                    }
                }
            }
        }
        
        // Get the pressure.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                    j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                    j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    p[idx_primitive_var] = V[d_num_species + d_dim.getValue()][idx_primitive_var];
                }
            }
        }
        
        // Compute the mass fractions.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_primitive_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                (ghostcell_dim_1_primitive_var + 1);
                        
                        Y[si][idx_primitive_var] = V[si][idx_primitive_var]/rho[idx_primitive_var];
                    }
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            data_mass_fractions,
            1);
        
        // Set the partial densities.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_conservative_var;
                        j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                        j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                (ghostcell_dim_1_conservative_var + 1);
                        
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                (ghostcell_dim_1_primitive_var + 1);
                        
                        Q[si][idx_conservative_var] = V[si][idx_primitive_var];
                    }
                }
            }
        }
        
        // Set the momentum.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                    j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                    j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    Q[d_num_species][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species][idx_primitive_var];
                    Q[d_num_species + 1][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species + 1][idx_primitive_var];
                    Q[d_num_species + 2][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species + 2][idx_primitive_var];
                }
            }
        }
        
        // Set the total energy.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                    j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                    j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    Q[d_num_species + d_dim.getValue()][idx_conservative_var] = rho[idx_primitive_var]*
                        (epsilon[idx_primitive_var] + double(1)/double(2)*(
                        V[d_num_species][idx_primitive_var]*V[d_num_species][idx_primitive_var] + 
                        V[d_num_species + 1][idx_primitive_var]*V[d_num_species + 1][idx_primitive_var] +
                        V[d_num_species + 2][idx_primitive_var]*V[d_num_species + 2][idx_primitive_var]));
                }
            }
        }
        
        /*
         * Convert primitive variables to conservative variables in the z-direction.
         */
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                Q[count_eqn] = conservative_variables[vi]->getPointer(2, di);
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                V[count_eqn] = primitive_variables[vi]->getPointer(2, di);
                count_eqn++;
            }
        }
        
        rho     = data_density->getPointer(2, 0);
        epsilon = data_internal_energy->getPointer(2, 0);
        p       = data_pressure->getPointer(2, 0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(1, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_primitive_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                ghostcell_dim_1_primitive_var;
                        
                        rho[idx_primitive_var] += V[si][idx_primitive_var];
                    }
                }
            }
        }
        
        // Get the pressure.
        for (int k = -num_ghosts_2_primitive_var;
                k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
                k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    p[idx_primitive_var] = V[d_num_species + d_dim.getValue()][idx_primitive_var];
                }
            }
        }
        
        // Compute the mass fractions.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_primitive_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                ghostcell_dim_1_primitive_var;
                        
                        Y[si][idx_primitive_var] = V[si][idx_primitive_var]/rho[idx_primitive_var];
                    }
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            data_mass_fractions,
            2);
        
        // Set the partial densities.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_conservative_var;
                    k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
                    k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                ghostcell_dim_1_conservative_var;
                        
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                ghostcell_dim_1_primitive_var;
                        
                        Q[si][idx_conservative_var] = V[si][idx_primitive_var];
                    }
                }
            }
        }
        
        // Set the momentum.
        for (int k = -num_ghosts_2_conservative_var;
                k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
                k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    Q[d_num_species][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species][idx_primitive_var];
                    Q[d_num_species + 1][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species + 1][idx_primitive_var];
                    Q[d_num_species + 2][idx_conservative_var] = rho[idx_primitive_var]*V[d_num_species + 2][idx_primitive_var];
                }
            }
        }
        
        // Set the total energy.
        for (int k = -num_ghosts_2_conservative_var;
                k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
                k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    Q[d_num_species + d_dim.getValue()][idx_conservative_var] = rho[idx_primitive_var]*
                        (epsilon[idx_primitive_var] + double(1)/double(2)*(
                        V[d_num_species][idx_primitive_var]*V[d_num_species][idx_primitive_var] + 
                        V[d_num_species + 1][idx_primitive_var]*V[d_num_species + 1][idx_primitive_var] +
                        V[d_num_species + 2][idx_primitive_var]*V[d_num_species + 2][idx_primitive_var]));
                }
            }
        }
    }
}


/*
 * Check whether the given cell conservative variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::checkCellDataOfConservativeVariablesBounded(
    HAMERS_SHARED_PTR<pdat::CellData<int> >& bounded_flag,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables)
{
    NULL_USE(bounded_flag);
    NULL_USE(conservative_variables);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelBasicUtilitiesFourEqnConservative::"
        << "checkCellDataOfConservativeVariablesBounded()\n"
        << "Method checkCellDataOfConservativeVariablesBounded()"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Check whether the given side conservative variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::checkSideDataOfConservativeVariablesBounded(
    HAMERS_SHARED_PTR<pdat::SideData<int> >& bounded_flag,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    
    /*
     * Get the dimensions of box that covers the interior of patch.
     */
    
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_flag = bounded_flag->getGhostCellWidth();
    
    const hier::IntVector num_ghosts_conservative_var = conservative_variables[0]->
        getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of conservative variables.
     */
    
    const hier::IntVector ghostcell_dims_conservative_var = conservative_variables[0]->
        getGhostBox().numberCells();
    
    /*
     * Check the size of variables.
     */
    
    if (static_cast<int>(conservative_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "checkSideDataOfConservativeVariablesBounded()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        const hier::IntVector interior_dims_conservative_var =
            conservative_variables[ei]->getBox().numberCells();
        
        if (interior_dims_conservative_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "checkSideDataOfConservativeVariablesBounded()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
                << std::endl);
        }
    }
    const hier::IntVector interior_dims_flag = bounded_flag->getBox().numberCells();
    if (interior_dims_flag != interior_dims)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "checkSideDataOfConservativeVariablesBounded()\n"
            << "The interior dimension of the flag does not match that of patch."
            << std::endl);
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_conservative_var != conservative_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "checkSideDataOfConservativeVariablesBounded()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_flag != num_ghosts_conservative_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "checkSideDataOfConservativeVariablesBounded()\n"
            << "The ghost cell width of the flag does not match that of conservative variables."
            << std::endl);
    }
    
    bounded_flag->fillAll(1);
    
    // Create the side data of density.
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_density(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_var));
    
    data_density->fillAll(double(0));
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    int* are_bounded = nullptr;
    
    std::vector<double*> Q;
    Q.resize(d_num_eqn);
    
    double* rho = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        
        /*
         * Check if conservative variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        rho = data_density->getPointer(0);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_conservative_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = i + num_ghosts_0_conservative_var;
                
                rho[idx_face] += Q[si][idx_face];
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_conservative_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = i + num_ghosts_0_conservative_var;
                
                const double Y = Q[si][idx_face]/rho[idx_face];
                
                if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
        // Check if density and total energy are bounded.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_conservative_var;
             i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_conservative_var;
            
            if (rho[idx_face] > double(0))
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
            
            if (Q[d_num_species + d_dim.getValue()][idx_face] > double(0))
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[1];
        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[0];
        
        /*
         * Check if conservative variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        rho = data_density->getPointer(0);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_conservative_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                    
                    rho[idx_face] += Q[si][idx_face];
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_conservative_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                    
                    const double Y = Q[si][idx_face]/rho[idx_face];
                    
                    if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        // Check if density and total energy are bounded.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_conservative_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                if (rho[idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (Q[d_num_species + d_dim.getValue()][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
        /*
         * Check if conservative variables in the y-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(1);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(1);
        }
        
        rho = data_density->getPointer(1);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                    
                    rho[idx_face] += Q[si][idx_face];
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                    
                    const double Y = Q[si][idx_face]/rho[idx_face];
                    
                    if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        // Check if density and total energy are bounded.
        for (int j = -num_ghosts_1_conservative_var;
             j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
             j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                if (rho[idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (Q[d_num_species + d_dim.getValue()][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        const int num_ghosts_1_conservative_var = num_ghosts_conservative_var[1];
        const int num_ghosts_2_conservative_var = num_ghosts_conservative_var[2];
        const int ghostcell_dim_0_conservative_var = ghostcell_dims_conservative_var[0];
        const int ghostcell_dim_1_conservative_var = ghostcell_dims_conservative_var[1];
        
        /*
         * Check if conservative variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        rho = data_density->getPointer(0);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_conservative_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                            (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                                ghostcell_dim_1_conservative_var;
                        
                        rho[idx_face] += Q[si][idx_face];
                    }
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_conservative_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                            (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                                ghostcell_dim_1_conservative_var;
                        
                        const double Y = Q[si][idx_face]/rho[idx_face];
                        
                        if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                        {
                            are_bounded[idx_face] &= 1;
                        }
                        else
                        {
                            are_bounded[idx_face] &= 0;
                        }
                    }
                }
            }
        }
        
        // Check if density and total energy are bounded.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_conservative_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    if (rho[idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_species + d_dim.getValue()][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        /*
         * Check if conservative variables in the y-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(1);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(1);
        }
        
        rho = data_density->getPointer(1);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_conservative_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                (ghostcell_dim_1_conservative_var + 1);
                        
                        rho[idx_face] += Q[si][idx_face];
                    }
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_conservative_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                (ghostcell_dim_1_conservative_var + 1);
                        
                        const double Y = Q[si][idx_face]/rho[idx_face];
                        
                        if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                        {
                            are_bounded[idx_face] &= 1;
                        }
                        else
                        {
                            are_bounded[idx_face] &= 0;
                        }
                    }
                }
            }
        }
        
        // Check if density and total energy are bounded.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    if (rho[idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_species + d_dim.getValue()][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        /*
         * Check if conservative variables in the z-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(2);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(2);
        }
        
        rho = data_density->getPointer(2);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_conservative_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                ghostcell_dim_1_conservative_var;
                        
                        rho[idx_face] += Q[si][idx_face];
                    }
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_conservative_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                            (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                                ghostcell_dim_1_conservative_var;
                        
                        const double Y = Q[si][idx_face]/rho[idx_face];
                        
                        if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                        {
                            are_bounded[idx_face] &= 1;
                        }
                        else
                        {
                            are_bounded[idx_face] &= 0;
                        }
                    }
                }
            }
        }
        
        // Check if density and total energy are bounded.
        for (int k = -num_ghosts_2_conservative_var;
             k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    if (rho[idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_species + d_dim.getValue()][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
    }
}


/*
 * Check whether the given cell primitive variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::checkCellDataOfPrimitiveVariablesBounded(
    HAMERS_SHARED_PTR<pdat::CellData<int> >& bounded_flag,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& primitive_variables)
{
    NULL_USE(bounded_flag);
    NULL_USE(primitive_variables);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelBasicUtilitiesFourEqnConservative::"
        << "checkCellDataOfPrimitiveVariablesBounded()\n"
        << "Method checkCellDataOfPrimitiveVariablesBounded()"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Check whether the given side primitive variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::checkSideDataOfPrimitiveVariablesBounded(
    HAMERS_SHARED_PTR<pdat::SideData<int> >& bounded_flag,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    
    /*
     * Get the dimensions of box that covers the interior of patch.
     */
    
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_flag = bounded_flag->getGhostCellWidth();
    
    const hier::IntVector num_ghosts_primitive_var = primitive_variables[0]->
        getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of primitive variables.
     */
    
    const hier::IntVector ghostcell_dims_primitive_var = primitive_variables[0]->
        getGhostBox().numberCells();
    
    /*
     * Check the size of variables.
     */
    
    if (static_cast<int>(primitive_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "checkSideDataOfPrimitiveVariablesBounded()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        const hier::IntVector interior_dims_primitive_var =
            primitive_variables[ei]->getBox().numberCells();
        
        if (interior_dims_primitive_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "checkSideDataOfPrimitiveVariablesBounded()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    const hier::IntVector interior_dims_flag = bounded_flag->getBox().numberCells();
    if (interior_dims_flag != interior_dims)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "checkSideDataOfPrimitiveVariablesBounded()\n"
            << "The interior dimension of the flag does not match that of patch."
            << std::endl);
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_primitive_var != primitive_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "checkSideDataOfPrimitiveVariablesBounded()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_flag != num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "checkSideDataOfPrimitiveVariablesBounded()\n"
            << "The ghost cell width of the flag does not match that of primitive variables."
            << std::endl);
    }
    
    bounded_flag->fillAll(1);
    
    // Create the side data of density.
    HAMERS_SHARED_PTR<pdat::SideData<double> > data_density(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_var));
    
    data_density->fillAll(double(0));
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    int* are_bounded = nullptr;
    
    std::vector<double*> V;
    V.resize(d_num_eqn);
    
    double* rho = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        
        /*
         * Check if primitive variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        rho = data_density->getPointer(0);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = i + num_ghosts_0_primitive_var;
                
                rho[idx_face] += V[si][idx_face];
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = i + num_ghosts_0_primitive_var;
                
                const double Y = V[si][idx_face]/rho[idx_face];
                
                if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
        // Check if density and pressure are bounded.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
             i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_primitive_var;
            
            if (rho[idx_face] > double(0))
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
            
            if (V[d_num_species + d_dim.getValue()][idx_face] > double(0))
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        const int num_ghosts_1_primitive_var = num_ghosts_primitive_var[1];
        const int ghostcell_dim_0_primitive_var = ghostcell_dims_primitive_var[0];
        
        /*
         * Check if primitive variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        rho = data_density->getPointer(0);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_primitive_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                    
                    rho[idx_face] += V[si][idx_face];
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_primitive_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                    
                    const double Y = V[si][idx_face]/rho[idx_face];
                    
                    if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        // Check if density and pressure are bounded.
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_primitive_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                if (rho[idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (V[d_num_species + d_dim.getValue()][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
        /*
         * Check if primitive variables in the y-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(1);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(1);
        }
        
        rho = data_density->getPointer(1);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                    
                    rho[idx_face] += V[si][idx_face];
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                    
                    const double Y = V[si][idx_face]/rho[idx_face];
                    
                    if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        // Check if density and pressure are bounded.
        for (int j = -num_ghosts_1_primitive_var;
             j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
             j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                if (rho[idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (V[d_num_species + d_dim.getValue()][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        const int num_ghosts_1_primitive_var = num_ghosts_primitive_var[1];
        const int num_ghosts_2_primitive_var = num_ghosts_primitive_var[2];
        const int ghostcell_dim_0_primitive_var = ghostcell_dims_primitive_var[0];
        const int ghostcell_dim_1_primitive_var = ghostcell_dims_primitive_var[1];
        
        /*
         * Check if primitive variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        rho = data_density->getPointer(0);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_primitive_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                            (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                                ghostcell_dim_1_primitive_var;
                        
                        rho[idx_face] += V[si][idx_face];
                    }
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_primitive_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                            (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                                ghostcell_dim_1_primitive_var;
                        
                        const double Y = V[si][idx_face]/rho[idx_face];
                        
                        if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                        {
                            are_bounded[idx_face] &= 1;
                        }
                        else
                        {
                            are_bounded[idx_face] &= 0;
                        }
                    }
                }
            }
        }
        
        // Check if density and pressure are bounded.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_primitive_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    if (rho[idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_species + d_dim.getValue()][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        /*
         * Check if primitive variables in the y-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(1);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(1);
        }
        
        rho = data_density->getPointer(1);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_primitive_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                (ghostcell_dim_1_primitive_var + 1);
                        
                        rho[idx_face] += V[si][idx_face];
                    }
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_primitive_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                (ghostcell_dim_1_primitive_var + 1);
                        
                        const double Y = V[si][idx_face]/rho[idx_face];
                        
                        if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                        {
                            are_bounded[idx_face] &= 1;
                        }
                        else
                        {
                            are_bounded[idx_face] &= 0;
                        }
                    }
                }
            }
        }
        
        // Check if density and pressure are bounded.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    if (rho[idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_species + d_dim.getValue()][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
        
        /*
         * Check if primitive variables in the z-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(2);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(2);
        }
        
        rho = data_density->getPointer(2);
        
        // Compute density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_primitive_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                ghostcell_dim_1_primitive_var;
                        
                        rho[idx_face] += V[si][idx_face];
                    }
                }
            }
        }
        
        // Check if mass fractions are bounded.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_primitive_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                ghostcell_dim_1_primitive_var;
                        
                        const double Y = V[si][idx_face]/rho[idx_face];
                        
                        if (Y > d_Y_bound_lo && Y < d_Y_bound_up)
                        {
                            are_bounded[idx_face] &= 1;
                        }
                        else
                        {
                            are_bounded[idx_face] &= 0;
                        }
                    }
                }
            }
        }
        
        // Check if density and pressure are bounded.
        for (int k = -num_ghosts_2_primitive_var;
             k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    if (rho[idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_species + d_dim.getValue()][idx_face] > double(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                }
            }
        }
    }
}


/*
 * Register the required derived variables for transformation between conservative
 * variables and characteristic variables.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables(
    const hier::IntVector& num_subghosts,
    const AVERAGING::TYPE& averaging_type)
{
    NULL_USE(num_subghosts);
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    // Check whether a patch is already registered.
    if (!flow_model_tmp->hasRegisteredPatch())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    d_proj_var_conservative_averaging_type = averaging_type;
}


/*
 * Register the required derived variables for transformation between primitive variables
 * and characteristic variables.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
    const hier::IntVector& num_subghosts,
    const AVERAGING::TYPE& averaging_type)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
    // Check whether a patch is already registered.
    if (!flow_model_tmp->hasRegisteredPatch())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    d_proj_var_primitive_averaging_type = averaging_type;
    
    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("DENSITY", num_subghosts));
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("SOUND_SPEED", num_subghosts));
    
    flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
}


/*
 * Get the number of projection variables for transformation between conservative
 * variables and characteristic variables.
 */
int
FlowModelBasicUtilitiesFourEqnConservative::getNumberOfProjectionVariablesForConservativeVariables() const
{
    TBOX_ERROR(d_object_name
        << ": FlowModelBasicUtilitiesFourEqnConservative::"
        << "getNumberOfProjectionVariablesForConservativeVariables()\n"
        << "Method getNumberOfProjectionVariablesForConservativeVariables()"
        << " is not yet implemented."
        << std::endl);
    
    return 0;
}


/*
 * Get the number of projection variables for transformation between primitive variables
 * and characteristic variables.
 */
int
FlowModelBasicUtilitiesFourEqnConservative::getNumberOfProjectionVariablesForPrimitiveVariables() const
{
    return d_num_species + 2;
}


/*
 * Compute the side data of the projection variables for transformation between conservative variables and
 * characteristic variables.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::computeSideDataOfProjectionVariablesForConservativeVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& projection_variables)
{
    NULL_USE(projection_variables);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelBasicUtilitiesFourEqnConservative::"
        << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
        << "Method computeSideDataOfProjectionVariablesForConservativeVariables()"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute the side data of the projection variables for transformation between primitive variables and
 * characteristic variables.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::computeSideDataOfProjectionVariablesForPrimitiveVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& projection_variables)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const hier::IntVector& num_ghosts = flow_model_tmp->getNumberOfGhostCells();
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    
    /*
     * Get the dimensions of the interior and ghost boxes.
     */
    
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    hier::Box ghost_box = interior_box;
    ghost_box.grow(num_ghosts);
    const hier::IntVector ghostcell_dims = ghost_box.numberCells();
    
    /*
     * Get the number of ghost cells and ghost cell dimension of projection variables.
     */
    
    const hier::IntVector num_ghosts_projection_var = projection_variables[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_projection_var =
        projection_variables[0]->getGhostBox().numberCells();
    
    /*
     * Check the size of variables.
     */
    
    if (static_cast<int>(projection_variables.size()) != d_num_species + 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
            << "There should be number of species projection plus two variables."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int vi = 0; vi < static_cast<int>(projection_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_projection_var =
            projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int vi = 1; vi < static_cast<int>(projection_variables.size()); vi++)
    {
        if (num_ghosts_projection_var != projection_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                << "The projection variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var > num_ghosts)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables have ghost cell width larger than that of conservative variables."
            << std::endl);
    }
    
    // Get the cell data of the variable partial densities.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_partial_densities =
        flow_model_tmp->getCellData("PARTIAL_DENSITIES");
    
    // Get the cell data of total density and sound speed.
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_density =
        flow_model_tmp->getCellData("DENSITY");
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_sound_speed =
        flow_model_tmp->getCellData("SOUND_SPEED");
    
    // Get the numbers of ghost cells and ghost cell dimensions of of total density and sound speed.
    const hier::IntVector& num_subghosts_density = data_density->getGhostCellWidth();
    const hier::IntVector subghostcell_dims_density = data_density->getGhostBox().numberCells();
    
    const hier::IntVector& num_subghosts_sound_speed = data_sound_speed->getGhostCellWidth();
    const hier::IntVector subghostcell_dims_sound_speed = data_sound_speed->getGhostBox().numberCells();
    
    if (num_ghosts_projection_var > num_subghosts_density)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables have ghost cell width larger than that of density."
            << std::endl);
    }
    
    if (num_ghosts_projection_var > num_subghosts_sound_speed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables have ghost cell width larger than that of sound speed."
            << std::endl);
    }
    
    // Get the pointers to the cell data of partial densities, total density and sound speed.
    std::vector<double*> rho_Y;
    rho_Y.reserve(d_num_species);
    for (int si = 0; si < d_num_species; si++)
    {
        rho_Y.push_back(data_partial_densities->getPointer(si));
    }
    double* rho = data_density->getPointer(0);
    double* c = data_sound_speed->getPointer(0);
    
    /*
     * Declare pointers to different data.
     */
    
    std::vector<double*> rho_Y_average;
    rho_Y_average.resize(d_num_species);
    double* rho_average = nullptr;
    double* c_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_0_projection_var = num_ghosts_projection_var[0];
        const int num_subghosts_0_density = num_subghosts_density[0];
        const int num_subghosts_0_sound_speed = num_subghosts_sound_speed[0];
        
        switch (d_proj_var_primitive_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_average[si] = projection_variables[si]->getPointer(0);
                }
                rho_average = projection_variables[d_num_species]->getPointer(0);
                c_average = projection_variables[d_num_species + 1]->getPointer(0);
                
                for (int si = 0; si < d_num_species; si++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_projection_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i + num_ghosts_0_projection_var;
                        const int idx_L = i - 1 + num_ghosts_0;
                        const int idx_R = i + num_ghosts_0;
                        
                        rho_Y_average[si][idx_face_x] = double(1)/double(2)*(rho_Y[si][idx_L] + rho_Y[si][idx_R]);
                    }
                }
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_projection_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i + num_ghosts_0_projection_var;
                    const int idx_density_L = i - 1 + num_subghosts_0_density;
                    const int idx_density_R = i + num_subghosts_0_density;
                    const int idx_sound_speed_L = i - 1 + num_subghosts_0_sound_speed;
                    const int idx_sound_speed_R = i + num_subghosts_0_sound_speed;
                    
                    rho_average[idx_face_x] = double(1)/double(2)*(rho[idx_density_L] + rho[idx_density_R]);
                    c_average[idx_face_x] = double(1)/double(2)*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesFourEqnConservative::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesFourEqnConservative::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Unknown d_proj_var_primitive_averaging_type given."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        const int ghostcell_dim_0 = ghostcell_dims[0];
        
        const int num_ghosts_0_projection_var = num_ghosts_projection_var[0];
        const int num_ghosts_1_projection_var = num_ghosts_projection_var[1];
        const int ghostcell_dim_0_projection_var = ghostcell_dims_projection_var[0];
        
        const int num_subghosts_0_density = num_subghosts_density[0];
        const int num_subghosts_1_density= num_subghosts_density[1];
        const int subghostcell_dim_0_density = subghostcell_dims_density[0];
        
        const int num_subghosts_0_sound_speed = num_subghosts_sound_speed[0];
        const int num_subghosts_1_sound_speed = num_subghosts_sound_speed[1];
        const int subghostcell_dim_0_sound_speed = subghostcell_dims_sound_speed[0];
        
        switch (d_proj_var_primitive_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_average[si] = projection_variables[si]->getPointer(0);
                }
                rho_average = projection_variables[d_num_species]->getPointer(0);
                c_average = projection_variables[d_num_species + 1]->getPointer(0);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = -num_ghosts_0_projection_var;
                             i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*(ghostcell_dim_0_projection_var + 1);
                            
                            const int idx_L = (i - 1 + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const int idx_R = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            rho_Y_average[si][idx_face_x] = double(1)/double(2)*(rho_Y[si][idx_L] + rho_Y[si][idx_R]);
                        }
                    }
                }
                
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_projection_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = (i + num_ghosts_0_projection_var) +
                            (j + num_ghosts_1_projection_var)*(ghostcell_dim_0_projection_var + 1);
                        
                        const int idx_density_L = (i - 1 + num_subghosts_0_density) +
                            (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                        
                        const int idx_density_R = (i + num_subghosts_0_density) +
                            (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                        
                        const int idx_sound_speed_L = (i - 1 + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_sound_speed_R = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        rho_average[idx_face_x] = double(1)/double(2)*(rho[idx_density_L] + rho[idx_density_R]);
                        c_average[idx_face_x] = double(1)/double(2)*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                    }
                }
                
                /*
                 * Compute the projection variables in the y-direction.
                 */
                
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_average[si] = projection_variables[si]->getPointer(1);
                }
                rho_average = projection_variables[d_num_species]->getPointer(1);
                c_average = projection_variables[d_num_species + 1]->getPointer(1);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int j = -num_ghosts_1_projection_var;
                         j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                         j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var;
                            
                            const int idx_B = (i + num_ghosts_0) +
                                (j - 1 + num_ghosts_1)*ghostcell_dim_0;
                            
                            const int idx_T = (i + num_ghosts_0) +
                                (j + num_ghosts_1)*ghostcell_dim_0;
                            
                            rho_Y_average[si][idx_face_y] = double(1)/double(2)*(rho_Y[si][idx_B] + rho_Y[si][idx_T]);
                        }
                    }
                }
                
                for (int j = -num_ghosts_1_projection_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = (i + num_ghosts_0_projection_var) +
                            (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var;
                        
                        const int idx_density_B = (i + num_subghosts_0_density) +
                            (j - 1 + num_subghosts_1_density)*subghostcell_dim_0_density;
                        
                        const int idx_density_T = (i + num_subghosts_0_density) +
                            (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                        
                        const int idx_sound_speed_B = (i + num_subghosts_0_sound_speed) +
                            (j - 1 + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_sound_speed_T = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        rho_average[idx_face_y] = double(1)/double(2)*(rho[idx_density_B] + rho[idx_density_T]);
                        c_average[idx_face_y] = double(1)/double(2)*(c[idx_sound_speed_B] + c[idx_sound_speed_T]);
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesFourEqnConservative::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesFourEqnConservative::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Unknown d_proj_var_primitive_averaging_type given."
                    << std::endl);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        const int num_ghosts_2 = num_ghosts[2];
        const int ghostcell_dim_0 = ghostcell_dims[0];
        const int ghostcell_dim_1 = ghostcell_dims[1];
        
        const int num_ghosts_0_projection_var = num_ghosts_projection_var[0];
        const int num_ghosts_1_projection_var = num_ghosts_projection_var[1];
        const int num_ghosts_2_projection_var = num_ghosts_projection_var[2];
        const int ghostcell_dim_0_projection_var = ghostcell_dims_projection_var[0];
        const int ghostcell_dim_1_projection_var = ghostcell_dims_projection_var[1];
        
        const int num_subghosts_0_density = num_subghosts_density[0];
        const int num_subghosts_1_density = num_subghosts_density[1];
        const int num_subghosts_2_density = num_subghosts_density[2];
        const int subghostcell_dim_0_density = subghostcell_dims_density[0];
        const int subghostcell_dim_1_density = subghostcell_dims_density[1];
        
        const int num_subghosts_0_sound_speed = num_subghosts_sound_speed[0];
        const int num_subghosts_1_sound_speed = num_subghosts_sound_speed[1];
        const int num_subghosts_2_sound_speed = num_subghosts_sound_speed[2];
        const int subghostcell_dim_0_sound_speed = subghostcell_dims_sound_speed[0];
        const int subghostcell_dim_1_sound_speed = subghostcell_dims_sound_speed[1];
        
        switch (d_proj_var_primitive_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_average[si] = projection_variables[si]->getPointer(0);
                }
                rho_average = projection_variables[d_num_species]->getPointer(0);
                c_average = projection_variables[d_num_species + 1]->getPointer(0);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = -num_ghosts_0_projection_var;
                                 i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                                 i++)
                            {
                                // Compute the linear indices.
                                const int idx_face_x = (i + num_ghosts_0_projection_var) +
                                    (j + num_ghosts_1_projection_var)*(ghostcell_dim_0_projection_var + 1) +
                                    (k + num_ghosts_2_projection_var)*(ghostcell_dim_0_projection_var + 1)*
                                        ghostcell_dim_1_projection_var;
                                
                                const int idx_L = (i - 1 + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                const int idx_R = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                rho_Y_average[si][idx_face_x] = double(1)/double(2)*(rho_Y[si][idx_L] + rho_Y[si][idx_R]);
                            }
                        }
                    }
                }
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = -num_ghosts_0_projection_var;
                             i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*(ghostcell_dim_0_projection_var + 1) +
                                (k + num_ghosts_2_projection_var)*(ghostcell_dim_0_projection_var + 1)*
                                    ghostcell_dim_1_projection_var;
                            
                            const int idx_density_L = (i - 1 + num_subghosts_0_density) +
                                (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_density_R = (i + num_subghosts_0_density) +
                                (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_sound_speed_L = (i - 1 + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_sound_speed_R = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            rho_average[idx_face_x] = double(1)/double(2)*(rho[idx_density_L] + rho[idx_density_R]);
                            c_average[idx_face_x] = double(1)/double(2)*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the y-direction.
                 */
                
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_average[si] = projection_variables[si]->getPointer(1);
                }
                rho_average = projection_variables[d_num_species]->getPointer(1);
                c_average = projection_variables[d_num_species + 1]->getPointer(1);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = -num_ghosts_1_projection_var;
                             j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                             j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_face_y = (i + num_ghosts_0_projection_var) +
                                    (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var +
                                    (k + num_ghosts_2_projection_var)*ghostcell_dim_0_projection_var*
                                        (ghostcell_dim_1_projection_var + 1);
                                
                                const int idx_B = (i + num_ghosts_0) +
                                    (j - 1 + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                const int idx_T = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                rho_Y_average[si][idx_face_y] = double(1)/double(2)*(rho_Y[si][idx_B] + rho_Y[si][idx_T]);
                            }
                        }
                    }
                }
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = -num_ghosts_1_projection_var;
                         j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                         j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var +
                                (k + num_ghosts_2_projection_var)*ghostcell_dim_0_projection_var*
                                    (ghostcell_dim_1_projection_var + 1);
                            
                            const int idx_density_B = (i + num_subghosts_0_density) +
                                (j - 1 + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_density_T = (i + num_subghosts_0_density) +
                                (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_sound_speed_B = (i + num_subghosts_0_sound_speed) +
                                (j - 1 + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_sound_speed_T = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            rho_average[idx_face_y] = double(1)/double(2)*(rho[idx_density_B] + rho[idx_density_T]);
                            c_average[idx_face_y] = double(1)/double(2)*(c[idx_sound_speed_B] + c[idx_sound_speed_T]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the z-direction.
                 */
                
                for (int si = 0; si < d_num_species; si++)
                {
                    rho_Y_average[si] = projection_variables[si]->getPointer(2);
                }
                rho_average = projection_variables[d_num_species]->getPointer(2);
                c_average = projection_variables[d_num_species + 1]->getPointer(2);
                
                for (int si = 0; si < d_num_species; si++)
                {
                    for (int k = -num_ghosts_2_projection_var;
                         k < interior_dim_2 + 1 + num_ghosts_2_projection_var;
                         k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx_face_z = (i + num_ghosts_0_projection_var) +
                                    (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var +
                                    (k + num_ghosts_2_projection_var)*ghostcell_dim_0_projection_var*
                                        ghostcell_dim_1_projection_var;
                                
                                const int idx_B = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k - 1 + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                const int idx_F = (i + num_ghosts_0) +
                                    (j + num_ghosts_1)*ghostcell_dim_0 +
                                    (k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                rho_Y_average[si][idx_face_z] = double(1)/double(2)*(rho_Y[si][idx_B] + rho_Y[si][idx_F]);
                            }
                        }
                    }
                }
                
                for (int k = -num_ghosts_2_projection_var;
                     k < interior_dim_2 + 1 + num_ghosts_2_projection_var;
                     k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var +
                                (k + num_ghosts_2_projection_var)*ghostcell_dim_0_projection_var*
                                    ghostcell_dim_1_projection_var;
                            
                            const int idx_density_B = (i + num_subghosts_0_density) +
                                (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k - 1 + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_density_F = (i + num_subghosts_0_density) +
                                (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                                (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                                    subghostcell_dim_1_density;
                            
                            const int idx_sound_speed_B = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k - 1 + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_sound_speed_F = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            rho_average[idx_face_z] = double(1)/double(2)*(rho[idx_density_B] + rho[idx_density_F]);
                            c_average[idx_face_z] = double(1)/double(2)*(c[idx_sound_speed_B] + c[idx_sound_speed_F]);
                        }
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesFourEqnConservative::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesFourEqnConservative::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Unknown d_proj_var_primitive_averaging_type given."
                    << std::endl);
            }
        }
    }
}


/*
 * Compute the side data of characteristic variables from conservative variables.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::computeSideDataOfCharacteristicVariablesFromConservativeVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& projection_variables,
    const int& idx_offset)
{
    NULL_USE(characteristic_variables);
    NULL_USE(conservative_variables);
    NULL_USE(projection_variables);
    NULL_USE(idx_offset);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelBasicUtilitiesFourEqnConservative::"
        << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
        << "Method computeSideDataOfCharacteristicVariablesFromConservativeVariables()"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute the side data of characteristic variables from primitive variables.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& primitive_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& projection_variables,
    const int& idx_offset)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    
    /*
     * Get the dimensions of box that covers the interior of patch.
     */
    
    const hier::IntVector interior_dims = patch.getBox().numberCells();
    
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_characteristic_var = characteristic_variables[0]->
        getGhostCellWidth();
    
    std::vector<hier::IntVector> num_ghosts_primitive_var;
    num_ghosts_primitive_var.reserve(static_cast<int>(primitive_variables.size()));
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        num_ghosts_primitive_var.push_back(primitive_variables[vi]->
            getGhostCellWidth());
    }
    
    const hier::IntVector num_ghosts_projection_var = projection_variables[0]->getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of characteristic and primitive variables.
     */
    
    const hier::IntVector ghostcell_dims_characteristic_var = characteristic_variables[0]->
        getGhostBox().numberCells();
    
    std::vector<hier::IntVector> ghostcell_dims_primitive_var;
    ghostcell_dims_primitive_var.reserve(static_cast<int>(primitive_variables.size()));
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        ghostcell_dims_primitive_var.push_back(primitive_variables[vi]->
            getGhostBox().numberCells());
    }
    
    /*
     * Check the size of variables.
     */
    
    if (static_cast<int>(characteristic_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(primitive_variables.size()) != 3)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    if (primitive_variables[0]->getDepth() != d_num_species ||
        primitive_variables[1]->getDepth() != d_dim.getValue() ||
        primitive_variables[2]->getDepth() != 1)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The depths of one or more primitive variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != d_num_species + 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "There should be number of species projection plus two variables."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        const hier::IntVector interior_dims_characteristic_var =
            characteristic_variables[ei]->getBox().numberCells();
        
        if (interior_dims_characteristic_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The interior dimension of the characteristic variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_primitive_var =
            primitive_variables[vi]->getBox().numberCells();
        
        if (interior_dims_primitive_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < d_num_species + 2; vi++)
    {
        const hier::IntVector interior_dims_projection_var = projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_characteristic_var != characteristic_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    for (int vi = 1; vi < d_num_species + 2; vi++)
    {
        if (num_ghosts_projection_var != projection_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The projection variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The ghost cell width of the projection variables does not match that of"
            << " characteristic variables."
            << std::endl);
    }
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        if (num_ghosts_primitive_var[vi] - num_ghosts_characteristic_var
                + (hier::IntVector::getOne(d_dim))*idx_offset < hier::IntVector::getZero(d_dim) ||
            num_ghosts_characteristic_var - num_ghosts_primitive_var[vi]
                + (hier::IntVector::getOne(d_dim))*(idx_offset + 1) > hier::IntVector::getZero(d_dim))
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The offset index is too large or the number of ghost of characteristic variable"
                << " is too large."
                << std::endl);
        }
    }
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<double*> V;
    V.reserve(d_num_eqn);
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        int depth = primitive_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            V.push_back(primitive_variables[vi]->getPointer(di));
        }
    }
    
    std::vector<double*> W;
    W.resize(d_num_eqn);
    
    std::vector<double*> rho_Y_average;
    rho_Y_average.resize(d_num_species);
    double* rho_average = nullptr;
    double* c_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_0_rho_Y = num_ghosts_primitive_var[0][0];
        const int num_ghosts_0_vel = num_ghosts_primitive_var[1][0];
        const int num_ghosts_0_p = num_ghosts_primitive_var[2][0];
        
        const int idx_offset_rho_Y = idx_offset;
        const int idx_offset_vel = idx_offset;
        const int idx_offset_p = idx_offset;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(0);
        }
        rho_average = projection_variables[d_num_species]->getPointer(0);
        c_average = projection_variables[d_num_species + 1]->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_characteristic_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                 i++)
            {
                // Compute the linear indices.
                const int idx_face = i + num_ghosts_0_characteristic_var;
                const int idx_rho_Y = i + idx_offset_rho_Y + num_ghosts_0_rho_Y;
                const int idx_p = i + idx_offset_p + num_ghosts_0_p;
                
                W[1 + si][idx_face] = V[si][idx_rho_Y] - rho_Y_average[si][idx_face]/
                    (rho_average[idx_face]*c_average[idx_face]*
                        c_average[idx_face])*V[d_num_species + 1][idx_p];
            }
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_characteristic_var;
             i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
             i++)
        {
            // Compute the linear indices.
            const int idx_face = i + num_ghosts_0_characteristic_var;
            const int idx_vel = i + idx_offset_vel + num_ghosts_0_vel;
            const int idx_p = i + idx_offset_p + num_ghosts_0_p;
            
            W[0][idx_face] = V[d_num_species][idx_vel] -
                double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 1][idx_p];
            
            W[d_num_species + 1][idx_face] = V[d_num_species][idx_vel] +
                double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 1][idx_p];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        
        const int num_ghosts_0_rho_Y = num_ghosts_primitive_var[0][0];
        const int num_ghosts_1_rho_Y = num_ghosts_primitive_var[0][1];
        const int ghostcell_dim_0_rho_Y = ghostcell_dims_primitive_var[0][0];
        
        const int num_ghosts_0_vel = num_ghosts_primitive_var[1][0];
        const int num_ghosts_1_vel = num_ghosts_primitive_var[1][1];
        const int ghostcell_dim_0_vel = ghostcell_dims_primitive_var[1][0];
        
        const int num_ghosts_0_p = num_ghosts_primitive_var[2][0];
        const int num_ghosts_1_p = num_ghosts_primitive_var[2][1];
        const int ghostcell_dim_0_p = ghostcell_dims_primitive_var[2][0];
        
        const int idx_offset_rho_Y = idx_offset;
        const int idx_offset_vel = idx_offset;
        const int idx_offset_p = idx_offset;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(0);
        }
        rho_average = projection_variables[d_num_species]->getPointer(0);
        c_average = projection_variables[d_num_species + 1]->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_characteristic_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1);
                    
                    const int idx_rho_Y = (i + idx_offset_rho_Y + num_ghosts_0_rho_Y) +
                        (j + num_ghosts_1_rho_Y)*ghostcell_dim_0_rho_Y;
                    
                    const int idx_p = (i + idx_offset_p + num_ghosts_0_p) +
                        (j + num_ghosts_1_p)*ghostcell_dim_0_p;
                    
                    W[1 + si][idx_face] = V[si][idx_rho_Y] - rho_Y_average[si][idx_face]/
                        (rho_average[idx_face]*c_average[idx_face]*c_average[idx_face])*
                            V[d_num_species + 2][idx_p];
                }
            }
        }
        
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_characteristic_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                 i++)
            {
                // Compute the linear indices.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1);
                
                const int idx_vel = (i + idx_offset_vel + num_ghosts_0_vel) +
                    (j + num_ghosts_1_vel)*ghostcell_dim_0_vel;
                
                const int idx_p = (i + idx_offset_p + num_ghosts_0_p) +
                    (j + num_ghosts_1_p)*ghostcell_dim_0_p;
                
                W[0][idx_face] = V[d_num_species][idx_vel] -
                    double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 2][idx_p];
                
                W[d_num_species + 1][idx_face] = V[d_num_species + 1][idx_vel];
                
                W[d_num_species + 2][idx_face] = V[d_num_species][idx_vel] +
                    double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 2][idx_p];
            }
        }
        
        /*
         * Compute the characteristic variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(1);
        }
        rho_average = projection_variables[d_num_species]->getPointer(1);
        c_average = projection_variables[d_num_species + 1]->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_characteristic_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var;
                    
                    const int idx_rho_Y = (i + num_ghosts_0_rho_Y) +
                        (j + idx_offset_rho_Y + num_ghosts_1_rho_Y)*ghostcell_dim_0_rho_Y;
                    
                    const int idx_p = (i + num_ghosts_0_p) +
                        (j + idx_offset_p + num_ghosts_1_p)*ghostcell_dim_0_p;
                    
                    W[1 + si][idx_face] = V[si][idx_rho_Y] - rho_Y_average[si][idx_face]/
                        (rho_average[idx_face]*c_average[idx_face]*c_average[idx_face])
                            *V[d_num_species + 2][idx_p];
                }
            }
        }
        
        for (int j = -num_ghosts_1_characteristic_var;
             j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
             j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var;
                
                const int idx_vel = (i + num_ghosts_0_vel) +
                    (j + idx_offset_vel + num_ghosts_1_vel)*ghostcell_dim_0_vel;
                
                const int idx_p = (i + num_ghosts_0_p) +
                    (j + idx_offset_p + num_ghosts_1_p)*ghostcell_dim_0_p;
                
                W[0][idx_face] = V[d_num_species + 1][idx_vel] -
                    double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 2][idx_p];
                
                W[d_num_species + 1][idx_face] = V[d_num_species][idx_vel];
                
                W[d_num_species + 2][idx_face] = V[d_num_species + 1][idx_vel] +
                    double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 2][idx_p];
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int num_ghosts_2_characteristic_var = num_ghosts_characteristic_var[2];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        const int ghostcell_dim_1_characteristic_var = ghostcell_dims_characteristic_var[1];
        
        const int num_ghosts_0_rho_Y = num_ghosts_primitive_var[0][0];
        const int num_ghosts_1_rho_Y = num_ghosts_primitive_var[0][1];
        const int num_ghosts_2_rho_Y = num_ghosts_primitive_var[0][2];
        const int ghostcell_dim_0_rho_Y = ghostcell_dims_primitive_var[0][0];
        const int ghostcell_dim_1_rho_Y = ghostcell_dims_primitive_var[0][1];
        
        const int num_ghosts_0_vel = num_ghosts_primitive_var[1][0];
        const int num_ghosts_1_vel = num_ghosts_primitive_var[1][1];
        const int num_ghosts_2_vel = num_ghosts_primitive_var[1][2];
        const int ghostcell_dim_0_vel = ghostcell_dims_primitive_var[1][0];
        const int ghostcell_dim_1_vel = ghostcell_dims_primitive_var[1][1];
        
        const int num_ghosts_0_p = num_ghosts_primitive_var[2][0];
        const int num_ghosts_1_p = num_ghosts_primitive_var[2][1];
        const int num_ghosts_2_p = num_ghosts_primitive_var[2][2];
        const int ghostcell_dim_0_p = ghostcell_dims_primitive_var[2][0];
        const int ghostcell_dim_1_p = ghostcell_dims_primitive_var[2][1];
        
        const int idx_offset_rho_Y = idx_offset;
        const int idx_offset_vel = idx_offset;
        const int idx_offset_p = idx_offset;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(0);
        }
        rho_average = projection_variables[d_num_species]->getPointer(0);
        c_average = projection_variables[d_num_species + 1]->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_characteristic_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_characteristic_var) +
                            (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1) +
                            (k + num_ghosts_2_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1)*
                                ghostcell_dim_1_characteristic_var;
                        
                        const int idx_rho_Y = (i + idx_offset_rho_Y + num_ghosts_0_rho_Y) +
                            (j + num_ghosts_1_rho_Y)*ghostcell_dim_0_rho_Y +
                            (k + num_ghosts_2_rho_Y)*ghostcell_dim_0_rho_Y*
                                ghostcell_dim_1_rho_Y;
                        
                        const int idx_p = (i + idx_offset_p + num_ghosts_0_p) +
                            (j + num_ghosts_1_p)*ghostcell_dim_0_p +
                            (k + num_ghosts_2_p)*ghostcell_dim_0_p*
                                ghostcell_dim_1_p;
                        
                        W[1 + si][idx_face] = V[si][idx_rho_Y] - rho_Y_average[si][idx_face]/
                            (rho_average[idx_face]*c_average[idx_face]*c_average[idx_face])*
                                V[d_num_species + 3][idx_p];
                    }
                }
            }
        }
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_characteristic_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1) +
                        (k + num_ghosts_2_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1)*
                            ghostcell_dim_1_characteristic_var;
                    
                    const int idx_vel = (i + idx_offset_vel + num_ghosts_0_vel) +
                        (j + num_ghosts_1_vel)*ghostcell_dim_0_vel +
                        (k + num_ghosts_2_vel)*ghostcell_dim_0_vel*
                            ghostcell_dim_1_vel;
                    
                    const int idx_p = (i + idx_offset_p + num_ghosts_0_p) +
                        (j + num_ghosts_1_p)*ghostcell_dim_0_p +
                        (k + num_ghosts_2_p)*ghostcell_dim_0_p*
                            ghostcell_dim_1_p;
                    
                    W[0][idx_face] = V[d_num_species][idx_vel] -
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 3][idx_p];
                    
                    W[d_num_species + 1][idx_face] = V[d_num_species + 1][idx_vel];
                    
                    W[d_num_species + 2][idx_face] = V[d_num_species + 2][idx_vel];
                    
                    W[d_num_species + 3][idx_face] = V[d_num_species][idx_vel] +
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 3][idx_p];
                }
            }
        }
        
        /*
         * Compute the characteristic variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(1);
        }
        rho_average = projection_variables[d_num_species]->getPointer(1);
        c_average = projection_variables[d_num_species + 1]->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_characteristic_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_characteristic_var) +
                            (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                            (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                                (ghostcell_dim_1_characteristic_var + 1);
                        
                        const int idx_rho_Y = (i + num_ghosts_0_rho_Y) +
                            (j + idx_offset_rho_Y + num_ghosts_1_rho_Y)*ghostcell_dim_0_rho_Y +
                            (k + num_ghosts_2_rho_Y)*ghostcell_dim_0_rho_Y*
                                ghostcell_dim_1_rho_Y;
                        
                        const int idx_p = (i + num_ghosts_0_p) +
                            (j + idx_offset_p + num_ghosts_1_p)*ghostcell_dim_0_p +
                            (k + num_ghosts_2_p)*ghostcell_dim_0_p*
                                ghostcell_dim_1_p;
                        
                        W[1 + si][idx_face] = V[si][idx_rho_Y] - rho_Y_average[si][idx_face]/
                            (rho_average[idx_face]*c_average[idx_face]*c_average[idx_face])*
                                V[d_num_species + 3][idx_p];
                    }
                }
            }
        }
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_characteristic_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            (ghostcell_dim_1_characteristic_var + 1);
                    
                    const int idx_vel = (i + num_ghosts_0_vel) +
                        (j + idx_offset_vel + num_ghosts_1_vel)*ghostcell_dim_0_vel +
                        (k + num_ghosts_2_vel)*ghostcell_dim_0_vel*
                            ghostcell_dim_1_vel;
                    
                    const int idx_p = (i + num_ghosts_0_p) +
                        (j + idx_offset_p + num_ghosts_1_p)*ghostcell_dim_0_p +
                        (k + num_ghosts_2_p)*ghostcell_dim_0_p*
                            ghostcell_dim_1_p;
                    
                    W[0][idx_face] = V[d_num_species + 1][idx_vel] -
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 3][idx_p];
                    
                    W[d_num_species + 1][idx_face] = V[d_num_species][idx_vel];
                    
                    W[d_num_species + 2][idx_face] = V[d_num_species + 2][idx_vel];
                    
                    W[d_num_species + 3][idx_face] = V[d_num_species + 1][idx_vel] +
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 3][idx_p];
                }
            }
        }
        
        /*
         * Compute the characteristic variables in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(2);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(2);
        }
        rho_average = projection_variables[d_num_species]->getPointer(2);
        c_average = projection_variables[d_num_species + 1]->getPointer(2);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_characteristic_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_characteristic_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face = (i + num_ghosts_0_characteristic_var) +
                            (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                            (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                                ghostcell_dim_1_characteristic_var;
                        
                        const int idx_rho_Y = (i + num_ghosts_0_rho_Y) +
                            (j + num_ghosts_1_rho_Y)*ghostcell_dim_0_rho_Y +
                            (k + idx_offset_rho_Y + num_ghosts_2_rho_Y)*ghostcell_dim_0_rho_Y*
                                ghostcell_dim_1_rho_Y;
                        
                        const int idx_p = (i + num_ghosts_0_p) +
                            (j + num_ghosts_1_p)*ghostcell_dim_0_p +
                            (k + idx_offset_p + num_ghosts_2_p)*ghostcell_dim_0_p*
                                ghostcell_dim_1_p;
                        
                        W[1 + si][idx_face] = V[si][idx_rho_Y] - rho_Y_average[si][idx_face]/
                            (rho_average[idx_face]*c_average[idx_face]*c_average[idx_face])*
                                V[d_num_species + 3][idx_p];
                    }
                }
            }
        }
        
        for (int k = -num_ghosts_2_characteristic_var;
             k < interior_dim_2 + 1 + num_ghosts_2_characteristic_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            ghostcell_dim_1_characteristic_var;
                    
                    const int idx_vel = (i + num_ghosts_0_vel) +
                        (j + num_ghosts_1_vel)*ghostcell_dim_0_vel +
                        (k + idx_offset_vel + num_ghosts_2_vel)*ghostcell_dim_0_vel*
                            ghostcell_dim_1_vel;
                    
                    const int idx_p = (i + num_ghosts_0_p) +
                        (j + num_ghosts_1_p)*ghostcell_dim_0_p +
                        (k + idx_offset_p + num_ghosts_2_p)*ghostcell_dim_0_p*
                            ghostcell_dim_1_p;
                    
                    W[0][idx_face] = V[d_num_species + 2][idx_vel] -
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 3][idx_p];
                    
                    W[d_num_species + 1][idx_face] = V[d_num_species][idx_vel];
                    
                    W[d_num_species + 2][idx_face] = V[d_num_species + 1][idx_vel];
                    
                    W[d_num_species + 3][idx_face] = V[d_num_species + 2][idx_vel] +
                        double(1)/(rho_average[idx_face]*c_average[idx_face])*V[d_num_species + 3][idx_p];
                }
            }
        }
    }
}


/*
 * Compute the side data of conservative variables from characteristic variables.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::computeSideDataOfConservativeVariablesFromCharacteristicVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& projection_variables)
{
    NULL_USE(conservative_variables);
    NULL_USE(characteristic_variables);
    NULL_USE(projection_variables);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelBasicUtilitiesFourEqnConservative::"
        << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
        << "Method computeSideDataOfConservativeVariablesFromCharacteristicVariables()"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Compute the side data of primitive variables from characteristic variables.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& projection_variables)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    const hier::Patch& patch = flow_model_tmp->getRegisteredPatch();
    
    /*
     * Get the dimensions of box that covers the interior of patch.
     */
    
    const hier::IntVector interior_dims = patch.getBox().numberCells();
    
    /*
     * Get the numbers of ghost cells of the variables.
     */
    
    const hier::IntVector num_ghosts_primitive_var = primitive_variables[0]->
        getGhostCellWidth();
    
    const hier::IntVector num_ghosts_characteristic_var = characteristic_variables[0]->
        getGhostCellWidth();
    
    const hier::IntVector num_ghosts_projection_var = projection_variables[0]->getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of characteristic variables.
     */
    
    const hier::IntVector ghostcell_dims_characteristic_var = characteristic_variables[0]->
        getGhostBox().numberCells();
    
    /*
     * Check the size of variables.
     */
    
    if (static_cast<int>(primitive_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(characteristic_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != d_num_species + 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "There should be number of species projection plus two variables."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        const hier::IntVector interior_dims_primitive_var =
            primitive_variables[ei]->getBox().numberCells();
        
        if (interior_dims_primitive_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        const hier::IntVector interior_dims_characteristic_var =
            characteristic_variables[ei]->getBox().numberCells();
        
        if (interior_dims_characteristic_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the characteristic variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < d_num_species + 2; vi++)
    {
        const hier::IntVector interior_dims_projection_var = projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_primitive_var != primitive_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_characteristic_var != characteristic_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    for (int vi = 1; vi < d_num_species + 2; vi++)
    {
        if (num_ghosts_projection_var != projection_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesFourEqnConservative::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The projection variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var != num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of"
            << " primitive variables."
            << std::endl);
    }
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesFourEqnConservative::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of"
            << " characteristic variables."
            << std::endl);
    }
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<double*> V;
    std::vector<double*> W;
    V.resize(d_num_eqn);
    W.resize(d_num_eqn);
    
    std::vector<double*> rho_Y_average;
    rho_Y_average.resize(d_num_species);
    double* rho_average = nullptr;
    double* c_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        
        /*
         * Compute the primitive variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(0);
        }
        rho_average = projection_variables[d_num_species]->getPointer(0);
        c_average = projection_variables[d_num_species + 1]->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_characteristic_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = i + num_ghosts_0_characteristic_var;
                
                V[si][idx_face] = -double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*W[0][idx_face] +
                    W[si + 1][idx_face] + double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*
                        W[d_num_eqn - 1][idx_face];
            }
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_characteristic_var;
             i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_characteristic_var;
            
            V[d_num_species][idx_face] = double(1)/double(2)*W[0][idx_face] +
                double(1)/double(2)*W[d_num_eqn - 1][idx_face];
            
            V[d_num_species + 1][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                W[0][idx_face] + double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                    W[d_num_eqn - 1][idx_face];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        
        /*
         * Compute the primitive variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(0);
        }
        rho_average = projection_variables[d_num_species]->getPointer(0);
        c_average = projection_variables[d_num_species + 1]->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_characteristic_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1);
                    
                    V[si][idx_face] = -double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*W[0][idx_face] +
                        W[si + 1][idx_face] + double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*
                            W[d_num_eqn - 1][idx_face];
                }
            }
        }
        
        for (int j = 0; j < interior_dim_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_ghosts_0_characteristic_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1);
                
                V[d_num_species][idx_face] = double(1)/double(2)*W[0][idx_face] +
                    double(1)/double(2)*W[d_num_eqn - 1][idx_face];
                
                V[d_num_species + 1][idx_face] = W[d_num_species + 1][idx_face];
                
                V[d_num_species + 2][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                    W[0][idx_face] + double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                        W[d_num_eqn - 1][idx_face];
            }
        }
        
        /*
         * Compute the primitive variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(1);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(1);
        }
        rho_average = projection_variables[d_num_species]->getPointer(1);
        c_average = projection_variables[d_num_species + 1]->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_characteristic_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var;
                    
                    V[si][idx_face] = -double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*W[0][idx_face] +
                        W[si + 1][idx_face] + double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*
                            W[d_num_eqn - 1][idx_face];
                }
            }
        }
        
        for (int j = -num_ghosts_1_characteristic_var;
             j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
             j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var;
                
                V[d_num_species][idx_face] = W[d_num_species + 1][idx_face];
                
                V[d_num_species + 1][idx_face] = double(1)/double(2)*W[0][idx_face] +
                    double(1)/double(2)*W[d_num_eqn - 1][idx_face];
                
                V[d_num_species + 2][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                    W[0][idx_face] + double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                        W[d_num_eqn - 1][idx_face];
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int num_ghosts_2_characteristic_var = num_ghosts_characteristic_var[2];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        const int ghostcell_dim_1_characteristic_var = ghostcell_dims_characteristic_var[1];
        
        /*
         * Compute the primitive variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(0);
        }
        rho_average = projection_variables[d_num_species]->getPointer(0);
        c_average = projection_variables[d_num_species + 1]->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = -num_ghosts_0_characteristic_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_characteristic_var) +
                            (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1) +
                            (k + num_ghosts_2_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1)*
                                ghostcell_dim_1_characteristic_var;
                        
                        V[si][idx_face] = -double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*W[0][idx_face] +
                            W[si + 1][idx_face] + double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*
                                W[d_num_eqn - 1][idx_face];
                    }
                }
            }
        }
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_ghosts_0_characteristic_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1) +
                        (k + num_ghosts_2_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1)*
                            ghostcell_dim_1_characteristic_var;
                    
                    V[d_num_species][idx_face] = double(1)/double(2)*W[0][idx_face] +
                        double(1)/double(2)*W[d_num_eqn - 1][idx_face];
                    
                    V[d_num_species + 1][idx_face] = W[d_num_species + 1][idx_face];
                    
                    V[d_num_species + 2][idx_face] = W[d_num_species + 2][idx_face];
                    
                    V[d_num_species + 3][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                        W[0][idx_face] + double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                            W[d_num_eqn - 1][idx_face];
                }
            }
        }
        
        /*
         * Compute the primitive variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(1);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(1);
        }
        rho_average = projection_variables[d_num_species]->getPointer(1);
        c_average = projection_variables[d_num_species + 1]->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1_characteristic_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                     j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_characteristic_var) +
                            (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                            (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                                (ghostcell_dim_1_characteristic_var + 1);
                        
                        V[si][idx_face] = -double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*W[0][idx_face] +
                            W[si + 1][idx_face] + double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*
                                W[d_num_eqn - 1][idx_face];
                    }
                }
            }
        }
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_characteristic_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                 j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            (ghostcell_dim_1_characteristic_var + 1);
                    
                    V[d_num_species][idx_face] = W[d_num_species + 1][idx_face];
                    
                    V[d_num_species + 1][idx_face] = double(1)/double(2)*W[0][idx_face] +
                        double(1)/double(2)*W[d_num_eqn - 1][idx_face];
                    
                    V[d_num_species + 2][idx_face] = W[d_num_species + 2][idx_face];
                    
                    V[d_num_species + 3][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                        W[0][idx_face] + double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*
                            W[d_num_eqn - 1][idx_face];
                }
            }
        }
        
        /*
         * Compute the primitive variables in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(2);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(2);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            rho_Y_average[si] = projection_variables[si]->getPointer(2);
        }
        rho_average = projection_variables[d_num_species]->getPointer(2);
        c_average = projection_variables[d_num_species + 1]->getPointer(2);
        
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_characteristic_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_characteristic_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_face = (i + num_ghosts_0_characteristic_var) +
                            (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                            (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                                ghostcell_dim_1_characteristic_var;
                        
                        V[si][idx_face] = -double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*W[0][idx_face] +
                            W[si + 1][idx_face] + double(1)/double(2)*rho_Y_average[si][idx_face]/c_average[idx_face]*
                                W[d_num_eqn - 1][idx_face];
                    }
                }
            }
        }
        
        for (int k = -num_ghosts_2_characteristic_var;
             k < interior_dim_2 + 1 + num_ghosts_2_characteristic_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            ghostcell_dim_1_characteristic_var;
                    
                    V[d_num_species][idx_face] = W[d_num_species + 1][idx_face];
                    
                    V[d_num_species + 1][idx_face] = W[d_num_species + 2][idx_face];
                    
                    V[d_num_species + 2][idx_face] = double(1)/double(2)*W[0][idx_face] +
                        double(1)/double(2)*W[d_num_eqn - 1][idx_face];
                    
                    V[d_num_species + 3][idx_face] = -double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*W[0][idx_face] +
                        double(1)/double(2)*rho_average[idx_face]*c_average[idx_face]*W[d_num_eqn - 1][idx_face];
                }
            }
        }
    }
}
