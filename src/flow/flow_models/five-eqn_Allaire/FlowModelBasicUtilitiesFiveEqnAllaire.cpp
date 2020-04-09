#include "flow/flow_models/five-eqn_Allaire/FlowModelBasicUtilitiesFiveEqnAllaire.hpp"

/*
 * Convert conservative variables to primitive variables.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::convertConservativeVariablesToPrimitiveVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    const int num_eqn = d_flow_model_tmp->getNumberOfEquations();
    
    /*
     * Get the dimensions of box that covers the interior of patch.
     */
    
    const hier::Box interior_box = primitive_variables[0]->getBox();
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
    
    int num_eqn_primitive_var = 0;
    int num_eqn_conservative_var = 0;
    
    /*
     * Check the size of variables.
     */
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        num_eqn_primitive_var += primitive_variables[vi]->getDepth();
    }
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        num_eqn_conservative_var += conservative_variables[vi]->getDepth();
    }
    
    if (!(num_eqn_primitive_var == num_eqn || num_eqn_primitive_var == num_eqn + 1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFiveEqnAllaire::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    
    if (!(num_eqn_conservative_var == num_eqn || num_eqn_conservative_var == num_eqn + 1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFiveEqnAllaire::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int vi = 1; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_primitive_var =
            primitive_variables[vi]->getBox().numberCells();
        
        if (interior_dims_primitive_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFiveEqnAllaire::"
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
                << ": FlowModelFiveEqnAllaire::"
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
                << ": FlowModelFiveEqnAllaire::"
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
                << ": FlowModelFiveEqnAllaire::"
                << "convertConservativeVariablesToPrimitiveVariables()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_primitive_var > num_ghosts_conservative_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFiveEqnAllaire::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "The ghost cell width of primitive variables is larger than that of conservative variables."
            << std::endl);
    }
    
    /*
     * Declare the pointers to the primitive variables and conservative variables.
     */
    
    std::vector<double*> V;
    V.resize(num_eqn_primitive_var);
    
    std::vector<double*> Q;
    Q.resize(num_eqn_conservative_var);
    
    int count_eqn = 0;
    
    /*
     * Convert conservative variables to primitive variables.
     */
    
    // Create the temporary side data.
    boost::shared_ptr<pdat::SideData<double> > data_density(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_internal_energy(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_pressure(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_mass_fractions(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_conservative_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_volume_fractions(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_conservative_var));
    
    data_density->fillAll(double(0));
    
    double* rho     = nullptr;
    double* epsilon = nullptr;
    double* p       = nullptr;
    
    std::vector<double*> Y;
    Y.resize(d_num_species);
    
    std::vector<double*> Z;
    Z.resize(d_num_species - 1);
    
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
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z[si] = data_volume_fractions->getPointer(0, si);
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
        
        // Get the volume fractions.
        for (int si = 0; si < d_num_species - 1; si++)
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
                
                Z[si][idx_conservative_var] = Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var];
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            data_mass_fractions,
            data_volume_fractions,
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
        
        // Set the volume fractions.
        if (num_eqn_primitive_var == num_eqn + 1)
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
                
                V[num_eqn][idx_primitive_var] = double(1);
            }
            
            for (int si = 0; si < d_num_species - 1; si++)
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
                    
                    V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var] = Z[si][idx_conservative_var];
                    V[num_eqn][idx_primitive_var] -= Z[si][idx_conservative_var];
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species - 1; si++)
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
                    
                    V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var] = Z[si][idx_conservative_var];
                }
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
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z[si] = data_volume_fractions->getPointer(0, si);
        }
        
        // Compute the mixture density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                 j < interior_dim_1 + num_ghosts_1_conservative_var;
                 j++)
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
        
        // Get the volume fractions.
        for (int si = 0; si < d_num_species - 1; si++)
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
                    
                    Z[si][idx_conservative_var] = Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var];
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            data_mass_fractions,
            data_volume_fractions,
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
        
        // Set the volume fractions.
        if (num_eqn_primitive_var == num_eqn + 1)
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
                    
                    V[num_eqn][idx_primitive_var] = double(1);
                }
            }
            
            for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var] = Z[si][idx_conservative_var];
                        V[num_eqn][idx_primitive_var] -= Z[si][idx_conservative_var];
                    }
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var] = Z[si][idx_conservative_var];
                    }
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
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z[si] = data_volume_fractions->getPointer(1, si);
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
        
        // Get the volume fractions.
        for (int si = 0; si < d_num_species - 1; si++)
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
                    
                    Z[si][idx_conservative_var] = Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var];
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            data_mass_fractions,
            data_volume_fractions,
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
        
        // Set the volume fractions.
        if (num_eqn_primitive_var == num_eqn + 1)
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
                    
                    V[num_eqn][idx_primitive_var] = double(1);
                }
            }
            
            for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var] = Z[si][idx_conservative_var];
                        V[num_eqn][idx_primitive_var] -= Z[si][idx_conservative_var];
                    }
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var] = Z[si][idx_conservative_var];
                    }
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
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z[si] = data_volume_fractions->getPointer(0, si);
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
        
        // Get the volume fractions.
        for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Z[si][idx_conservative_var] = Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var];
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
            data_volume_fractions,
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
        
        // Set the volume fractions.
        if (num_eqn_primitive_var == num_eqn + 1)
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
                        
                        V[num_eqn][idx_primitive_var] = double(1);
                    }
                }
            }
            
            for (int si = 0; si < d_num_species - 1; si++)
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
                            
                            V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var] = Z[si][idx_conservative_var];
                            V[num_eqn][idx_primitive_var] -= Z[si][idx_conservative_var];
                        }
                    }
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species - 1; si++)
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
                            
                            V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var] = Z[si][idx_conservative_var];
                        }
                    }
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
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z[si] = data_volume_fractions->getPointer(1, si);
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
        
        // Get the volume fractions.
        for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Z[si][idx_conservative_var] = Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var];
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
            data_volume_fractions,
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
        
        // Set the volume fractions.
        if (num_eqn_primitive_var == num_eqn + 1)
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
                        
                        V[num_eqn][idx_primitive_var] = double(1);
                    }
                }
            }
            
            for (int si = 0; si < d_num_species - 1; si++)
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
                            
                            V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var] = Z[si][idx_conservative_var];
                            V[num_eqn][idx_primitive_var] -= Z[si][idx_conservative_var];
                        }
                    }
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species - 1; si++)
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
                            
                            V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var] = Z[si][idx_conservative_var];
                        }
                    }
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
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z[si] = data_volume_fractions->getPointer(2, si);
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
        
        // Get the volume fractions.
        for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Z[si][idx_conservative_var] = Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var];
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
            data_volume_fractions,
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
        
        // Set the volume fractions.
        if (num_eqn_primitive_var == num_eqn + 1)
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
                        
                        V[num_eqn][idx_primitive_var] = double(1);
                    }
                }
            }
            
            for (int si = 0; si < d_num_species - 1; si++)
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
                            
                            V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var] = Z[si][idx_conservative_var];
                            V[num_eqn][idx_primitive_var] -= Z[si][idx_conservative_var];
                        }
                    }
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species - 1; si++)
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
                            
                            V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var] = Z[si][idx_conservative_var];
                        }
                    }
                }
            }
        }
    }
}


/*
 * Convert primitive variables to conservative variables.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::convertPrimitiveVariablesToConservativeVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    const int num_eqn = d_flow_model_tmp->getNumberOfEquations();
    
    /*
     * Get the dimensions of box that covers the interior of patch.
     */
    
    const hier::Box interior_box = conservative_variables[0]->getBox();
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
    
    int num_eqn_conservative_var = 0;
    int num_eqn_primitive_var = 0;
    
    /*
     * Check the size of variables.
     */
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        num_eqn_conservative_var += conservative_variables[vi]->getDepth();
    }
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        num_eqn_primitive_var += primitive_variables[vi]->getDepth();
    }
    
    if (!(num_eqn_conservative_var == num_eqn || num_eqn_conservative_var == num_eqn + 1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFiveEqnAllaire::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    
    if (!(num_eqn_primitive_var == num_eqn || num_eqn_primitive_var == num_eqn + 1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFiveEqnAllaire::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int vi = 1; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        const hier::IntVector interior_dims_conservative_var =
            conservative_variables[vi]->getBox().numberCells();
        
        if (interior_dims_conservative_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFiveEqnAllaire::"
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
                << ": FlowModelFiveEqnAllaire::"
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
                << ": FlowModelFiveEqnAllaire::"
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
                << ": FlowModelFiveEqnAllaire::"
                << "convertPrimitiveVariablesToConservativeVariables()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_conservative_var > num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFiveEqnAllaire::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "The ghost cell width of conservative variables is larger than that of primitive variables."
            << std::endl);
    }
    
    /*
     * Declare the pointers to the conservative variables and primitive variables.
     */
    
    std::vector<double*> Q;
    Q.resize(num_eqn_conservative_var);
    
    std::vector<double*> V;
    V.resize(num_eqn_primitive_var);
    
    int count_eqn = 0;
    
    /*
     * Convert primitive variables to conservative variables.
     */
    
    // Create the temporary side data.
    boost::shared_ptr<pdat::SideData<double> > data_density(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_pressure(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_internal_energy(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_mass_fractions(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_volume_fractions(
        new pdat::SideData<double>(interior_box, d_num_species - 1, num_ghosts_primitive_var));
    
    data_density->fillAll(double(0));
    
    double* rho     = nullptr;
    double* p       = nullptr;
    double* epsilon = nullptr;
    
    std::vector<double*> Y;
    Y.resize(d_num_species);
    
    std::vector<double*> Z;
    Z.resize(d_num_species - 1);
    
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
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z[si] = data_volume_fractions->getPointer(0, si);
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
        
        // Get the volume fractions.
        for (int si = 0; si < d_num_species - 1; si++)
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
                
                Z[si][idx_primitive_var] = V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var];
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
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
        
        // Set the volume fractions.
        if (num_eqn_conservative_var == num_eqn + 1)
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
                
                Q[num_eqn][idx_conservative_var] = double(1);
            }
            
            for (int si = 0; si < d_num_species - 1; si++)
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
                    
                    Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var] = Z[si][idx_primitive_var];
                    Q[num_eqn][idx_conservative_var] -= Z[si][idx_primitive_var];
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species - 1; si++)
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
                    
                    Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var] = Z[si][idx_primitive_var];
                }
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
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z[si] = data_volume_fractions->getPointer(0, si);
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
        
        // Get the volume fractions.
        for (int si = 0; si < d_num_species - 1; si++)
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
                    
                    Z[si][idx_primitive_var] = V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var];
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
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
        
        // Set the volume fractions.
        if (num_eqn_conservative_var == num_eqn + 1)
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
                    
                    Q[num_eqn][idx_conservative_var] = double(1);
                }
            }
            
            for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var] = Z[si][idx_primitive_var];
                        Q[num_eqn][idx_conservative_var] -= Z[si][idx_primitive_var];
                    }
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var] = Z[si][idx_primitive_var];
                    }
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
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z[si] = data_volume_fractions->getPointer(1, si);
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
        
        // Get the volume fractions.
        for (int si = 0; si < d_num_species - 1; si++)
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
                    
                    Z[si][idx_primitive_var] = V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var];
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
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
        
        // Set the volume fractions.
        if (num_eqn_conservative_var == num_eqn + 1)
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
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                    
                    Q[num_eqn][idx_conservative_var] = double(1);
                }
            }
            
            for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var] = Z[si][idx_primitive_var];
                        Q[num_eqn][idx_conservative_var] -= Z[si][idx_primitive_var];
                    }
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var] = Z[si][idx_primitive_var];
                    }
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
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z[si] = data_volume_fractions->getPointer(0, si);
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
        
        // Get the volume fractions.
        for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Z[si][idx_primitive_var] = V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var];
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
            data_volume_fractions,
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
        
        // Set the volume fractions.
        if (num_eqn_conservative_var == num_eqn + 1)
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
                        
                        Q[num_eqn][idx_conservative_var] = double(1);
                    }
                }
            }
            
            for (int si = 0; si < d_num_species - 1; si++)
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
                            
                            Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var] = Z[si][idx_primitive_var];
                            Q[num_eqn][idx_conservative_var] -= Z[si][idx_primitive_var];
                        }
                    }
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species - 1; si++)
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
                            
                            Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var] = Z[si][idx_primitive_var];
                        }
                    }
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
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z[si] = data_volume_fractions->getPointer(1, si);
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
        
        // Get the volume fractions.
        for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Z[si][idx_primitive_var] = V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var];
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
            data_volume_fractions,
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
        
        // Set the volume fractions.
        if (num_eqn_conservative_var == num_eqn + 1)
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
                        
                        Q[num_eqn][idx_conservative_var] = double(1);
                    }
                }
            }
            
            for (int si = 0; si < d_num_species - 1; si++)
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
                            
                            Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var] = Z[si][idx_primitive_var];
                            Q[num_eqn][idx_conservative_var] -= Z[si][idx_primitive_var];
                        }
                    }
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species - 1; si++)
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
                            
                            Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var] = Z[si][idx_primitive_var];
                        }
                    }
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
        
        for (int si = 0; si < d_num_species - 1; si++)
        {
            Z[si] = data_volume_fractions->getPointer(1, si);
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
        
        // Get the volume fractions.
        for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Z[si][idx_primitive_var] = V[d_num_species + d_dim.getValue() + 1 + si][idx_primitive_var];
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
            data_volume_fractions,
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
        
        // Set the volume fractions.
        if (num_eqn_conservative_var == num_eqn + 1)
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
                        
                        Q[num_eqn][idx_conservative_var] = double(1);
                    }
                }
            }
            
            for (int si = 0; si < d_num_species - 1; si++)
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
                            
                            Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var] = Z[si][idx_primitive_var];
                            Q[num_eqn][idx_conservative_var] -= Z[si][idx_primitive_var];
                        }
                    }
                }
            }
        }
        else
        {
            for (int si = 0; si < d_num_species - 1; si++)
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
                            
                            Q[d_num_species + d_dim.getValue() + 1 + si][idx_conservative_var] = Z[si][idx_primitive_var];
                        }
                    }
                }
            }
        }
    }
}


/*
 * Check whether the given cell conservative variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::checkCellDataOfConservativeVariablesBounded(
    boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables)
{
    // NEED IMPLEMENTATION!
}


/*
 * Check whether the given side conservative variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::checkSideDataOfConservativeVariablesBounded(
    boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    const int num_eqn = d_flow_model_tmp->getNumberOfEquations();
    
    /*
     * Get the dimensions of box that covers the interior of patch.
     */
    
    const hier::Box interior_box = bounded_flag->getBox();
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
    
    if (!(static_cast<int>(conservative_variables.size()) == num_eqn ||
          static_cast<int>(conservative_variables.size()) - 1 == num_eqn))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFiveEqnAllaire::"
            << "checkSideDataOfConservativeVariablesBounded()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        const hier::IntVector interior_dims_conservative_var =
            conservative_variables[ei]->getBox().numberCells();
        
        if (interior_dims_conservative_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFiveEqnAllaire::"
                << "checkSideDataOfConservativeVariablesBounded()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < num_eqn; ei++)
    {
        if (num_ghosts_conservative_var != conservative_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFiveEqnAllaire::"
                << "checkSideDataOfConservativeVariablesBounded()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_flag != num_ghosts_conservative_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFiveEqnAllaire::"
            << "checkSideDataOfConservativeVariablesBounded()\n"
            << "The ghost cell width of the flag does not match that of conservative variables."
            << std::endl);
    }
    
    bounded_flag->fillAll(1);
    
    // Create the side data for last volume fraction.
    boost::shared_ptr<pdat::SideData<double> > data_last_volume_fraction(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_var));
    
    data_last_volume_fraction->fillAll(double(1));
    
    // Create the side data of density.
    boost::shared_ptr<pdat::SideData<double> > data_density(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_conservative_var));
    
    data_density->fillAll(double(0));
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    int* are_bounded = nullptr;
    
    std::vector<double*> Q;
    Q.resize(num_eqn);
    
    double* Z_last = nullptr;
    
    double* rho = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_conservative_var = num_ghosts_conservative_var[0];
        
        /*
         * Check if conservative variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        Z_last = data_last_volume_fraction->getPointer(0);
        
        rho = data_density->getPointer(0);
        
        // Compute last volume fraction and check if volume fractions are bounded.
        for (int si = 0; si < d_num_species - 1; si++)
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
                
                Z_last[idx_face] -= Q[d_num_species + d_dim.getValue() + 1 + si][idx_face];
                
                if (Q[d_num_species + d_dim.getValue() + 1 + si][idx_face] > d_Z_bound_lo &&
                    Q[d_num_species + d_dim.getValue() + 1 + si][idx_face] < d_Z_bound_up)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
        // Check if last volume fraction is bounded.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_conservative_var;
             i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_conservative_var;
            
            if (Z_last[idx_face] > d_Z_bound_lo && Z_last[idx_face] < d_Z_bound_up)
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
        }
        
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
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        Z_last = data_last_volume_fraction->getPointer(0);
        
        rho = data_density->getPointer(0);
        
        // Compute last volume fraction and check if volume fractions are bounded.
        for (int si = 0; si < d_num_species - 1; si++)
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
                    
                    Z_last[idx_face] -= Q[d_num_species + d_dim.getValue() + 1 + si][idx_face];
                    
                    if (Q[d_num_species + d_dim.getValue() + 1 + si][idx_face] > d_Z_bound_lo &&
                        Q[d_num_species + d_dim.getValue() + 1 + si][idx_face] < d_Z_bound_up)
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
        
        // Check if last volume fraction is bounded.
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
                
                if (Z_last[idx_face] > d_Z_bound_lo && Z_last[idx_face] < d_Z_bound_up)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
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
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(1);
        }
        
        Z_last = data_last_volume_fraction->getPointer(1);
        
        rho = data_density->getPointer(1);
        
        // Compute last volume fraction and check if volume fractions are bounded.
        for (int si = 0; si < d_num_species - 1; si++)
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
                    
                    Z_last[idx_face] -= Q[d_num_species + d_dim.getValue() + 1 + si][idx_face];
                    
                    if (Q[d_num_species + d_dim.getValue() + 1 + si][idx_face] > d_Z_bound_lo &&
                        Q[d_num_species + d_dim.getValue() + 1 + si][idx_face] < d_Z_bound_up)
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
        
        // Check if last volume fraction is bounded.
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
                
                if (Z_last[idx_face] > d_Z_bound_lo && Z_last[idx_face] < d_Z_bound_up)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
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
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        Z_last = data_last_volume_fraction->getPointer(0);
        
        rho = data_density->getPointer(0);
        
        // Compute last volume fraction and check if volume fractions are bounded.
        for (int si = 0; si < d_num_species - 1; si++)
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
                        const int idx_face = (i + num_ghosts_0_conservative_var) +
                            (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                            (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                                ghostcell_dim_1_conservative_var;
                        
                        Z_last[idx_face] -= Q[d_num_species + d_dim.getValue() + 1 + si][idx_face];
                        
                        if (Q[d_num_species + d_dim.getValue() + 1 + si][idx_face] > d_Z_bound_lo &&
                            Q[d_num_species + d_dim.getValue() + 1 + si][idx_face] < d_Z_bound_up)
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
        
        // Check if last volume fraction is bounded.
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
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    if (Z_last[idx_face] > d_Z_bound_lo && Z_last[idx_face] < d_Z_bound_up)
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
                        // Compute the linear index.
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
                        // Compute the linear index.
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
                    // Compute the linear index.
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
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(1);
        }
        
        Z_last = data_last_volume_fraction->getPointer(1);
        
        rho = data_density->getPointer(1);
        
        // Compute last volume fraction and check if volume fractions are bounded.
        for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Z_last[idx_face] -= Q[d_num_species + d_dim.getValue() + 1 + si][idx_face];
                        
                        if (Q[d_num_species + d_dim.getValue() + 1 + si][idx_face] > d_Z_bound_lo &&
                            Q[d_num_species + d_dim.getValue() + 1 + si][idx_face] < d_Z_bound_up)
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
        
        // Check if last volume fraction is bounded.
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
                    
                    if (Z_last[idx_face] > d_Z_bound_lo && Z_last[idx_face] < d_Z_bound_up)
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
                        // Compute the linear index.
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
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(2);
        }
        
        Z_last = data_last_volume_fraction->getPointer(2);
        
        rho = data_density->getPointer(2);
        
        // Compute last volume fraction and check if volume fractions are bounded.
        for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Z_last[idx_face] -= Q[d_num_species + d_dim.getValue() + 1 + si][idx_face];
                        
                        if (Q[d_num_species + d_dim.getValue() + 1 + si][idx_face] > d_Z_bound_lo &&
                            Q[d_num_species + d_dim.getValue() + 1 + si][idx_face] < d_Z_bound_up)
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
        
        // Check if last volume fraction is bounded.
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
                    
                    if (Z_last[idx_face] > d_Z_bound_lo && Z_last[idx_face] < d_Z_bound_up)
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
                        // Compute the linear index.
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
FlowModelBasicUtilitiesFiveEqnAllaire::checkCellDataOfPrimitiveVariablesBounded(
    boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables)
{
    // NEED IMPLEMENTATION!
}


/*
 * Check whether the given side primitive variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::checkSideDataOfPrimitiveVariablesBounded(
    boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables)
{
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    const int num_eqn = d_flow_model_tmp->getNumberOfEquations();
    
    /*
     * Get the dimensions of box that covers the interior of patch.
     */
    
    const hier::Box interior_box = bounded_flag->getBox();
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
    
    if (!(static_cast<int>(primitive_variables.size()) == num_eqn ||
          static_cast<int>(primitive_variables.size()) - 1 == num_eqn))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFiveEqnAllaire::"
            << "checkSideDataOfPrimitiveVariablesBounded()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        const hier::IntVector interior_dims_primitive_var =
            primitive_variables[ei]->getBox().numberCells();
        
        if (interior_dims_primitive_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFiveEqnAllaire::"
                << "checkSideDataOfPrimitiveVariablesBounded()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < num_eqn; ei++)
    {
        if (num_ghosts_primitive_var != primitive_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelFiveEqnAllaire::"
                << "checkSideDataOfPrimitiveVariablesBounded()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_flag != num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelFiveEqnAllaire::"
            << "checkSideDataOfPrimitiveVariablesBounded()\n"
            << "The ghost cell width of the flag does not match that of primitive variables."
            << std::endl);
    }
    
    bounded_flag->fillAll(1);
    
    // Create the side data for volume fractions.
    boost::shared_ptr<pdat::SideData<double> > data_volume_fractions(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_var));
    
    data_volume_fractions->fill(double(1), d_num_species - 1);
    
    // Create the side data of density.
    boost::shared_ptr<pdat::SideData<double> > data_density(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_var));
    
    data_density->fillAll(double(0));
    
    // Create the side data of pressure.
    boost::shared_ptr<pdat::SideData<double> > data_pressure(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_var));
    
    // Create the side data of mass fractions.
    boost::shared_ptr<pdat::SideData<double> > data_mass_fractions(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_var));
    
    // Create the side data of partial derivatives.
    boost::shared_ptr<pdat::SideData<double> > data_gruneisen_parameter(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_var));
    
    boost::shared_ptr<pdat::SideData<double> > data_partial_pressure_partial_partial_densities(
        new pdat::SideData<double>(interior_box, d_num_species, num_ghosts_primitive_var));
    
    // Create the side data of square of sound speed.
    boost::shared_ptr<pdat::SideData<double> > data_sound_speed_sq(
        new pdat::SideData<double>(interior_box, 1, num_ghosts_primitive_var));
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    int* are_bounded = nullptr;
    
    std::vector<double*> V;
    V.resize(num_eqn);
    
    std::vector<double*> Z;
    Z.resize(num_eqn);
    
    double* rho = nullptr;
    double* p = nullptr;
    
    std::vector<double*> Y;
    Y.resize(num_eqn);
    
    double* Gamma;
    
    std::vector<double*> Psi;
    Psi.resize(d_num_species);
    
    double* c_sq;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_primitive_var = num_ghosts_primitive_var[0];
        
        /*
         * Check if primitive variables in the x-direction are bounded.
         */
        
        are_bounded = bounded_flag->getPointer(0);
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            Z[si] = data_volume_fractions->getPointer(0, si);
        }
        
        rho = data_density->getPointer(0);
        p = data_pressure->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(0, si);
        }
        
        Gamma = data_gruneisen_parameter->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Psi[si] = data_partial_pressure_partial_partial_densities->getPointer(0, si);
        }
        
        c_sq = data_sound_speed_sq->getPointer(0);
        
        // Get the volume fractions, compute last volume fraction and check if volume fractions are bounded.
        for (int si = 0; si < d_num_species - 1; si++)
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
                
                Z[si][idx_face] = V[d_num_species + d_dim.getValue() + 1 + si][idx_face];
                Z[d_num_species - 1][idx_face] -= Z[si][idx_face];
                
                if (Z[si][idx_face] > d_Z_bound_lo &&
                    Z[si][idx_face] < d_Z_bound_up)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
        // Check if last volume fraction is bounded.
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
             i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_primitive_var;
            
            if (Z[d_num_species - 1][idx_face] > d_Z_bound_lo &&
                Z[d_num_species - 1][idx_face] < d_Z_bound_up)
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
        }
        
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
                
                Y[si][idx_face] = V[si][idx_face]/rho[idx_face];
                
                if (Y[si][idx_face] > d_Y_bound_lo &&
                    Y[si][idx_face] < d_Y_bound_up)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
        // Check if partial densities are bounded.
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
                
                if (V[si][idx_face] > double(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
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
            const int idx_face = i + num_ghosts_0_primitive_var;
            
            p[idx_face] = V[d_num_species + d_dim.getValue()][idx_face];
        }
        
        // Check if sound speed is real.
        
        d_equation_of_state_mixing_rules->computeGruneisenParameter(
            data_gruneisen_parameter,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
            0);
        
        d_equation_of_state_mixing_rules->computePressureDerivativeWithPartialDensities(
            data_partial_pressure_partial_partial_densities,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
            0);
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
             i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_primitive_var;
            
            c_sq[idx_face] = Gamma[idx_face]*p[idx_face]/rho[idx_face];
        }
        
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
                
                c_sq[idx_face] += Y[si][idx_face]*Psi[si][idx_face];
            }
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_ghosts_0_primitive_var;
             i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_primitive_var;
            
            if (c_sq[idx_face] > double(0))
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
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            Z[si] = data_volume_fractions->getPointer(0, si);
        }
        
        rho = data_density->getPointer(0);
        p = data_pressure->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(0, si);
        }
        
        Gamma = data_gruneisen_parameter->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Psi[si] = data_partial_pressure_partial_partial_densities->getPointer(0, si);
        }
        
        c_sq = data_sound_speed_sq->getPointer(0);
        
        // Get the volume fractions, compute last volume fraction and check if volume fractions are bounded.
        for (int si = 0; si < d_num_species - 1; si++)
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
                    
                    Z[si][idx_face] = V[d_num_species + d_dim.getValue() + 1 + si][idx_face];
                    Z[d_num_species - 1][idx_face] -= Z[si][idx_face];
                    
                    if (Z[si][idx_face] > d_Z_bound_lo &&
                        Z[si][idx_face] < d_Z_bound_up)
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
        
        // Check if last volume fraction is bounded.
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
                
                if (Z[d_num_species - 1][idx_face] > d_Z_bound_lo &&
                    Z[d_num_species - 1][idx_face] < d_Z_bound_up)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
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
                    
                    Y[si][idx_face] = V[si][idx_face]/rho[idx_face];
                    
                    if (Y[si][idx_face] > d_Y_bound_lo &&
                        Y[si][idx_face] < d_Y_bound_up)
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
        
        // Check if partial densities are bounded.
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
                    
                    if (V[si][idx_face] > double(0))
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
                const int idx_face = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                p[idx_face] = V[d_num_species + d_dim.getValue()][idx_face];
            }
        }
        
        // Check if sound speed is real.
        
        d_equation_of_state_mixing_rules->computeGruneisenParameter(
            data_gruneisen_parameter,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
            0);
        
        d_equation_of_state_mixing_rules->computePressureDerivativeWithPartialDensities(
            data_partial_pressure_partial_partial_densities,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
            0);
        
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
                
                c_sq[idx_face] = Gamma[idx_face]*p[idx_face]/rho[idx_face];
            }
        }
        
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
                        
                        c_sq[idx_face] += Y[si][idx_face]*Psi[si][idx_face];
                }
            }
        }
        
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
                
                if (c_sq[idx_face] > double(0))
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
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(1);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            Z[si] = data_volume_fractions->getPointer(1, si);
        }
        
        rho = data_density->getPointer(1);
        p = data_pressure->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(1, si);
        }
        
        Gamma = data_gruneisen_parameter->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Psi[si] = data_partial_pressure_partial_partial_densities->getPointer(1, si);
        }
        
        c_sq = data_sound_speed_sq->getPointer(1);
        
        // Get the volume fractions, compute last volume fraction and check if volume fractions are bounded.
        for (int si = 0; si < d_num_species - 1; si++)
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
                    
                    Z[si][idx_face] = V[d_num_species + d_dim.getValue() + 1 + si][idx_face];
                    Z[d_num_species - 1][idx_face] -= Z[si][idx_face];
                    
                    if (Z[si][idx_face] > d_Z_bound_lo &&
                        Z[si][idx_face] < d_Z_bound_up)
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
        
        // Check if last volume fraction is bounded.
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
                
                if (Z[d_num_species - 1][idx_face] > d_Z_bound_lo &&
                    Z[d_num_species - 1][idx_face] < d_Z_bound_up)
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
            }
        }
        
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
                    
                    Y[si][idx_face] = V[si][idx_face]/rho[idx_face];
                    
                    if (Y[si][idx_face] > d_Y_bound_lo &&
                        Y[si][idx_face] < d_Y_bound_up)
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
        
        // Check if partial densities are bounded.
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
                    
                    if (V[si][idx_face] > double(0))
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
                const int idx_face = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                p[idx_face] = V[d_num_species + d_dim.getValue()][idx_face];
            }
        }
        
        // Check if sound speed is real.
        
        d_equation_of_state_mixing_rules->computeGruneisenParameter(
            data_gruneisen_parameter,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
            1);
        
        d_equation_of_state_mixing_rules->computePressureDerivativeWithPartialDensities(
            data_partial_pressure_partial_partial_densities,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
            1);
        
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
                
                c_sq[idx_face] = Gamma[idx_face]*p[idx_face]/rho[idx_face];
            }
        }
        
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
                    
                    c_sq[idx_face] = Gamma[idx_face]*p[idx_face]/rho[idx_face];
                }
            }
        }
        
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
                
                if (c_sq[idx_face] > double(0))
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
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(0);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            Z[si] = data_volume_fractions->getPointer(0, si);
        }
        
        rho = data_density->getPointer(0);
        p = data_pressure->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(0, si);
        }
        
        Gamma = data_gruneisen_parameter->getPointer(0);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Psi[si] = data_partial_pressure_partial_partial_densities->getPointer(0, si);
        }
        
        c_sq = data_sound_speed_sq->getPointer(0);
        
        // Get the volume fractions, compute last volume fraction and check if volume fractions are bounded.
        for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Z[si][idx_face] = V[d_num_species + d_dim.getValue() + 1 + si][idx_face];
                        Z[d_num_species - 1][idx_face] -= Z[si][idx_face];
                        
                        if (Z[si][idx_face] > d_Z_bound_lo &&
                            Z[si][idx_face] < d_Z_bound_up)
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
        
        // Check if last volume fraction is bounded.
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
                    
                    if (Z[d_num_species - 1][idx_face] > d_Z_bound_lo &&
                        Z[d_num_species - 1][idx_face] < d_Z_bound_up)
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
                        
                        Y[si][idx_face] = V[si][idx_face]/rho[idx_face];
                        
                        if (Y[si][idx_face] > d_Y_bound_lo &&
                            Y[si][idx_face] < d_Y_bound_up)
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
        
        // Check if partial densities are bounded.
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
                        
                        if (V[si][idx_face] > double(0))
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
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    p[idx_face] = V[d_num_species + d_dim.getValue()][idx_face];
                }
            }
        }
        
        // Check if sound speed is real.
        
        d_equation_of_state_mixing_rules->computeGruneisenParameter(
            data_gruneisen_parameter,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
            0);
        
        d_equation_of_state_mixing_rules->computePressureDerivativeWithPartialDensities(
            data_partial_pressure_partial_partial_densities,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
            0);
        
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
                    
                    c_sq[idx_face] = Gamma[idx_face]*p[idx_face]/rho[idx_face];
                }
            }
        }
        
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
                        
                        c_sq[idx_face] += Y[si][idx_face]*Psi[si][idx_face];
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
                for (int i = -num_ghosts_0_primitive_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    if (c_sq[idx_face] > double(0))
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
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(1);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            Z[si] = data_volume_fractions->getPointer(1, si);
        }
        
        rho = data_density->getPointer(1);
        p = data_pressure->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(1, si);
        }
        
        Gamma = data_gruneisen_parameter->getPointer(1);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Psi[si] = data_partial_pressure_partial_partial_densities->getPointer(1, si);
        }
        
        c_sq = data_sound_speed_sq->getPointer(1);
        
        // Get the volume fractions, compute last volume fraction and check if volume fractions are bounded.
        for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Z[si][idx_face] = V[d_num_species + d_dim.getValue() + 1 + si][idx_face];
                        Z[d_num_species - 1][idx_face] -= Z[si][idx_face];
                        
                        if (Z[si][idx_face] > d_Z_bound_lo &&
                            Z[si][idx_face] < d_Z_bound_up)
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
        
        // Check if last volume fraction is bounded.
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
                    
                    if (Z[d_num_species - 1][idx_face] > d_Z_bound_lo &&
                        Z[d_num_species - 1][idx_face] < d_Z_bound_up)
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
                        
                        Y[si][idx_face] = V[si][idx_face]/rho[idx_face];
                        
                        if (Y[si][idx_face] > d_Y_bound_lo &&
                            Y[si][idx_face] < d_Y_bound_up)
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
        
        // Check if partial densities are bounded.
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
                        
                        if (V[si][idx_face] > double(0))
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
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    p[idx_face] = V[d_num_species + d_dim.getValue()][idx_face];
                }
            }
        }
        
        // Check if sound speed is real.
        
        d_equation_of_state_mixing_rules->computeGruneisenParameter(
            data_gruneisen_parameter,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
            1);
        
        d_equation_of_state_mixing_rules->computePressureDerivativeWithPartialDensities(
            data_partial_pressure_partial_partial_densities,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
            1);
        
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
                    
                    c_sq[idx_face] = Gamma[idx_face]*p[idx_face]/rho[idx_face];
                }
            }
        }
        
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
                        
                        c_sq[idx_face] = Gamma[idx_face]*p[idx_face]/rho[idx_face];
                    }
                }
            }
        }
        
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
                    
                    if (c_sq[idx_face] > double(0))
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
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
            V[ei] = primitive_variables[ei]->getPointer(2);
        }
        
        for (int si = 0; si < d_num_species; si++)
        {
            Z[si] = data_volume_fractions->getPointer(2, si);
        }
        
        rho = data_density->getPointer(2);
        p = data_pressure->getPointer(2);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Y[si] = data_mass_fractions->getPointer(2, si);
        }
        
        Gamma = data_gruneisen_parameter->getPointer(2);
        
        for (int si = 0; si < d_num_species; si++)
        {
            Psi[si] = data_partial_pressure_partial_partial_densities->getPointer(2, si);
        }
        
        c_sq = data_sound_speed_sq->getPointer(2);
        
        // Get the volume fractions, compute last volume fraction and check if volume fractions are bounded.
        for (int si = 0; si < d_num_species - 1; si++)
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
                        
                        Z[si][idx_face] = V[d_num_species + d_dim.getValue() + 1 + si][idx_face];
                        Z[d_num_species - 1][idx_face] -= Z[si][idx_face];
                        
                        if (Z[si][idx_face] > d_Z_bound_lo &&
                            Z[si][idx_face] < d_Z_bound_up)
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
        
        // Check if last volume fraction is bounded.
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
                    
                    if (Z[d_num_species - 1][idx_face] > d_Z_bound_lo &&
                        Z[d_num_species - 1][idx_face] < d_Z_bound_up)
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
                        
                        Y[si][idx_face] = V[si][idx_face]/rho[idx_face];
                        
                        if (Y[si][idx_face] > d_Y_bound_lo &&
                            Y[si][idx_face] < d_Y_bound_up)
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
        
        // Check if partial densities are bounded.
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
                        
                        if (V[si][idx_face] > double(0))
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
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    p[idx_face] = V[d_num_species + d_dim.getValue()][idx_face];
                }
            }
        }
        
        // Check if sound speed is real.
        
        d_equation_of_state_mixing_rules->computeGruneisenParameter(
            data_gruneisen_parameter,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
            2);
        
        d_equation_of_state_mixing_rules->computePressureDerivativeWithPartialDensities(
            data_partial_pressure_partial_partial_densities,
            data_density,
            data_pressure,
            data_mass_fractions,
            data_volume_fractions,
            2);
        
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
                    
                    c_sq[idx_face] = Gamma[idx_face]*p[idx_face]/rho[idx_face];
                }
            }
        }
        
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
                        
                        c_sq[idx_face] = Gamma[idx_face]*p[idx_face]/rho[idx_face];
                    }
                }
            }
        }
        
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
                    
                    if (c_sq[idx_face] > double(0))
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
 * Get the number of projection variables for transformation between conservative
 * variables and characteristic variables.
 */
int
FlowModelBasicUtilitiesFiveEqnAllaire::getNumberOfProjectionVariablesForConservativeVariables() const
{
    TBOX_ERROR(d_object_name
        << ": FlowModelBasicUtilitiesFiveEqnAllaire::"
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
FlowModelBasicUtilitiesFiveEqnAllaire::getNumberOfProjectionVariablesForPrimitiveVariables() const
{
    return d_num_species + 2;
}


/*
 * Compute the side data of the projection variables for transformation between conservative variables and
 * characteristic variables.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::computeSideDataOfProjectionVariablesForConservativeVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables)
{
}


/*
 * Compute the side data of the projection variables for transformation between primitive variables and
 * characteristic variables.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::computeSideDataOfProjectionVariablesForPrimitiveVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables)
{
}


/*
 * Compute the side data of characteristic variables from conservative variables.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::computeSideDataOfCharacteristicVariablesFromConservativeVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables,
    const int& idx_offset)
{
}


/*
 * Compute the side data of characteristic variables from primitive variables.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables,
    const int& idx_offset)
{
}


/*
 * Compute the side data of conservative variables from characteristic variables.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::computeSideDataOfConservativeVariablesFromCharacteristicVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables)
{
}


/*
 * Compute the side data of primitive variables from characteristic variables.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables)
{
}
