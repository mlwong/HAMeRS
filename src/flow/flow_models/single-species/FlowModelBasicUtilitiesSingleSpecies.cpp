#include "flow/flow_models/single-species/FlowModelBasicUtilitiesSingleSpecies.hpp"

/*
 * Convert conservative variables to primitive variables.
 */
void
FlowModelBasicUtilitiesSingleSpecies::convertConservativeVariablesToPrimitiveVariables(
    const std::vector<const Real*>& conservative_variables,
    const std::vector<Real*>& primitive_variables)
{
    const std::vector<const Real*>& Q = conservative_variables;
    const std::vector<Real*>&       V = primitive_variables;
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    if (static_cast<int>(Q.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "convertLocalCellDataPointersConservativeVariablesToPrimitiveVariables()\n"
            << "Number of elements in conservative variables is not correct."
            << std::endl);
    }
    
    if (static_cast<int>(V.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "convertLocalCellDataPointersConservativeVariablesToPrimitiveVariables()\n"
            << "Number of elements in primitive variables is not correct."
            << std::endl);
    }
#endif
    
    // Get the pointers to the momentum components.
    std::vector<const Real*> m_ptr;
    m_ptr.reserve(d_dim.getValue());
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        m_ptr.push_back(Q[1 + di]);
    }
    
    /*
     * Get the thermodynamic properties of the species.
     */
    
    const int num_thermo_properties = d_equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<Real> thermo_properties;
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    // Compute the specific internal energy.
    Real epsilon = Real(0);
    if (d_dim == tbox::Dimension(1))
    {
        epsilon = (*Q[2] - Real(1)/Real(2)*(*Q[1])*(*Q[1])/(*Q[0]))/(*Q[0]);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        epsilon = (*Q[3] - Real(1)/Real(2)*((*Q[1])*(*Q[1]) + (*Q[2])*(*Q[2]))/(*Q[0]))/(*Q[0]);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        epsilon = (*Q[4] - Real(1)/Real(2)*((*Q[1])*(*Q[1]) + (*Q[2])*(*Q[2]) + (*Q[3])*(*Q[3]))/(*Q[0]))/(*Q[0]);
    }
    
    // Compute the pressure.
    const Real p = d_equation_of_state_mixing_rules->getEquationOfState()->
        getPressure(
            Q[0],
            &epsilon,
            thermo_properties_const_ptr);
    
    // Convert the conservative variables to primitive variables.
    *V[0] = *Q[0];
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        *V[1 + di] = (*Q[1 + di])/(*Q[0]);
    }
    *V[1 + d_dim.getValue()] = p;
}


/*
 * Convert conservative variables to primitive variables.
 */
void
FlowModelBasicUtilitiesSingleSpecies::convertConservativeVariablesToPrimitiveVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables)
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
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    
    if (d_num_eqn_conservative_var != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "convertConservativeVariablesToPrimitiveVariables()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_primitive_var > num_ghosts_conservative_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "convertConservativeVariablesToPrimitiveVariables()\n"
            << "The ghost cell width of primitive variables is larger than that of conservative variables."
            << std::endl);
    }
    
    /*
     * Declare the pointers to the primitive variables and conservative variables.
     */
    
    std::vector<Real*> V;
    V.resize(d_num_eqn_primitive_var);
    
    std::vector<Real*> Q;
    Q.resize(d_num_eqn_conservative_var);
    
    int count_eqn = 0;
    
    /*
     * Convert conservative variables to primitive variables.
     */
    
    // Create the temporary side data.
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_density(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_internal_energy(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_pressure(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_conservative_var));
    
    Real* rho     = nullptr;
    Real* epsilon = nullptr;
    Real* p       = nullptr;
    
    /*
     * Get the thermodynamic properties of the species.
     */
    
    const int num_thermo_properties = d_equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<Real> thermo_properties;
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
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
        
        // Get the density.
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear index.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            rho[idx_conservative_var] = Q[0][idx_conservative_var];
        }
        
        // Compute the internal energy.
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear index.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            epsilon[idx_conservative_var] = (Q[1 + d_dim.getValue()][idx_conservative_var] -
                Real(1)/Real(2)*(Q[1][idx_conservative_var]*Q[1][idx_conservative_var])/
                rho[idx_conservative_var])/rho[idx_conservative_var];
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            thermo_properties_const_ptr,
            0);
        
        // Set the density.
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            V[0][idx_primitive_var] = Q[0][idx_conservative_var];
        }
        
        // Set the velocity.
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            V[1][idx_primitive_var] = Q[1][idx_conservative_var]/rho[idx_conservative_var];
        }
        
        // Set the pressure.
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            
            V[1 + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
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
        
        // Get the density.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_conservative_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                    i++)
            {
                // Compute the linear index.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                rho[idx_conservative_var] = Q[0][idx_conservative_var];
            }
        }
        
        // Compute the internal energy.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_conservative_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                    i++)
            {
                // Compute the linear index.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                epsilon[idx_conservative_var] = (Q[1 + d_dim.getValue()][idx_conservative_var] -
                    Real(1)/Real(2)*(Q[1][idx_conservative_var]*Q[1][idx_conservative_var] +
                    Q[2][idx_conservative_var]*Q[2][idx_conservative_var])/
                    rho[idx_conservative_var])/rho[idx_conservative_var];
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            thermo_properties_const_ptr,
            0);
        
        // Set the density.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_primitive_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                V[0][idx_primitive_var] = Q[0][idx_conservative_var];
            }
        }
        
        // Set the velocity.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_primitive_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                V[1][idx_primitive_var] = Q[1][idx_conservative_var]/rho[idx_conservative_var];
                V[2][idx_primitive_var] = Q[2][idx_conservative_var]/rho[idx_conservative_var];
            }
        }
        
        // Set the pressure.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_primitive_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                V[1 + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
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
        
        // Get the density.
        for (int j = -num_ghosts_1_conservative_var;
                j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                rho[idx_conservative_var] = Q[0][idx_conservative_var];
            }
        }
        
        // Compute the internal energy.
        for (int j = -num_ghosts_1_conservative_var;
                j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                epsilon[idx_conservative_var] = (Q[1 + d_dim.getValue()][idx_conservative_var] -
                    Real(1)/Real(2)*(Q[1][idx_conservative_var]*Q[1][idx_conservative_var] +
                    Q[2][idx_conservative_var]*Q[2][idx_conservative_var])/
                    rho[idx_conservative_var])/rho[idx_conservative_var];
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            thermo_properties_const_ptr,
            1);
        
        // Set the density.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                V[0][idx_primitive_var] = Q[0][idx_conservative_var];
            }
        }
        
        // Set the velocity.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                V[1][idx_primitive_var] = Q[1][idx_conservative_var]/rho[idx_conservative_var];
                V[2][idx_primitive_var] = Q[2][idx_conservative_var]/rho[idx_conservative_var];
            }
        }
        
        // Set the pressure.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                V[1 + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
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
        
        // Get the density.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_conservative_var;
                        i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                        i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    rho[idx_conservative_var] = Q[0][idx_conservative_var];
                }
            }
        }
        
        // Compute the internal energy.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_conservative_var;
                        i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                        i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    epsilon[idx_conservative_var] = (Q[1 + d_dim.getValue()][idx_conservative_var] -
                        Real(1)/Real(2)*(Q[1][idx_conservative_var]*Q[1][idx_conservative_var] +
                        Q[2][idx_conservative_var]*Q[2][idx_conservative_var] +
                        Q[3][idx_conservative_var]*Q[3][idx_conservative_var])/
                        rho[idx_conservative_var])/rho[idx_conservative_var];
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            thermo_properties_const_ptr,
            0);
        
        // Set the density.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    V[0][idx_primitive_var] = Q[0][idx_conservative_var];
                }
            }
        }
        
        // Set the velocity.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    V[1][idx_primitive_var] = Q[1][idx_conservative_var]/rho[idx_conservative_var];
                    V[2][idx_primitive_var] = Q[2][idx_conservative_var]/rho[idx_conservative_var];
                    V[3][idx_primitive_var] = Q[3][idx_conservative_var]/rho[idx_conservative_var];
                }
            }
        }
        
        // Set the pressure.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    V[1 + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
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
        
        // Get the density.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                    j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                    j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    rho[idx_conservative_var] = Q[0][idx_conservative_var];
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
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    epsilon[idx_conservative_var] = (Q[1 + d_dim.getValue()][idx_conservative_var] -
                        Real(1)/Real(2)*(Q[1][idx_conservative_var]*Q[1][idx_conservative_var] +
                        Q[2][idx_conservative_var]*Q[2][idx_conservative_var] +
                        Q[3][idx_conservative_var]*Q[3][idx_conservative_var])/
                        rho[idx_conservative_var])/rho[idx_conservative_var];
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            thermo_properties_const_ptr,
            1);
        
        // Set the density.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                    j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                    j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    V[0][idx_primitive_var] = Q[0][idx_conservative_var];
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
                HAMERS_PRAGMA_SIMD
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
                    
                    V[1][idx_primitive_var] = Q[1][idx_conservative_var]/rho[idx_conservative_var];
                    V[2][idx_primitive_var] = Q[2][idx_conservative_var]/rho[idx_conservative_var];
                    V[3][idx_primitive_var] = Q[3][idx_conservative_var]/rho[idx_conservative_var];
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
                HAMERS_PRAGMA_SIMD
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
                    
                    V[1 + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
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
        
        // Get the density.
        for (int k = -num_ghosts_2_conservative_var;
                k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
                k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    rho[idx_conservative_var] = Q[0][idx_conservative_var];
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
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    epsilon[idx_conservative_var] = (Q[1 + d_dim.getValue()][idx_conservative_var] -
                        Real(1)/Real(2)*(Q[1][idx_conservative_var]*Q[1][idx_conservative_var] +
                        Q[2][idx_conservative_var]*Q[2][idx_conservative_var] +
                        Q[3][idx_conservative_var]*Q[3][idx_conservative_var])/
                        rho[idx_conservative_var])/rho[idx_conservative_var];
                }
            }
        }
        
        // Compute the pressure.
        d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
            data_pressure,
            data_density,
            data_internal_energy,
            thermo_properties_const_ptr,
            2);
        
        // Set the density.
        for (int k = -num_ghosts_2_primitive_var;
                k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
                k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    V[0][idx_primitive_var] = Q[0][idx_conservative_var];
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
                HAMERS_PRAGMA_SIMD
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
                    
                    V[1][idx_primitive_var] = Q[1][idx_conservative_var]/rho[idx_conservative_var];
                    V[2][idx_primitive_var] = Q[2][idx_conservative_var]/rho[idx_conservative_var];
                    V[3][idx_primitive_var] = Q[3][idx_conservative_var]/rho[idx_conservative_var];
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
                HAMERS_PRAGMA_SIMD
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
                    
                    V[1 + d_dim.getValue()][idx_primitive_var] = p[idx_conservative_var];
                }
            }
        }
    }
}


/*
 * Convert primitive variables to conservative variables.
 */
void
FlowModelBasicUtilitiesSingleSpecies::convertPrimitiveVariablesToConservativeVariables(
    const std::vector<const Real*>& primitive_variables,
    const std::vector<Real*>& conservative_variables)
{
    const std::vector<const Real*>& V = primitive_variables;
    const std::vector<Real*>&       Q = conservative_variables;
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    if (static_cast<int>(V.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "Number of elements in primitive variables is not correct."
            << std::endl);
    }
    
    if (static_cast<int>(Q.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "Number of elements in conservative variables is not correct."
            << std::endl);
    }
#endif
    
    /*
     * Get the thermodynamic properties of the species.
     */
    
    const int num_thermo_properties = d_equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<Real> thermo_properties;
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    // Compute the total energy.
    const Real epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
        getInternalEnergy(
            V[0],
            V[1 + d_dim.getValue()],
            thermo_properties_const_ptr);
    
    Real E = Real(0);
    if (d_dim == tbox::Dimension(1))
    {
        E = (*V[0])*(epsilon + Real(1)/Real(2)*(*V[1])*(*V[1]));
    }
    else if (d_dim == tbox::Dimension(2))
    {
        E = (*V[0])*(epsilon + Real(1)/Real(2)*((*V[1])*(*V[1]) + (*V[2])*(*V[2])));
    }
    else if (d_dim == tbox::Dimension(3))
    {
        E = (*V[0])*(epsilon + Real(1)/Real(2)*((*V[1])*(*V[1]) + (*V[2])*(*V[2]) + (*V[3])*(*V[3])));
    }
    
    // Convert the primitive variables to conservative variables.
    *Q[0] = *V[0];
    for (int di = 0; di < d_dim.getValue(); di++)
    {
        *Q[1 + di] = (*V[0])*(*V[1 + di]);
    }
    *Q[1 + d_dim.getValue()] = E;
}


/*
 * Convert primitive variables to conservative variables.
 */
void
FlowModelBasicUtilitiesSingleSpecies::convertPrimitiveVariablesToConservativeVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables)
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
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    
    if (d_num_eqn_primitive_var != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "convertPrimitiveVariablesToConservativeVariables()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_conservative_var > num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "convertPrimitiveVariablesToConservativeVariables()\n"
            << "The ghost cell width of conservative variables is larger than that of primitive variables."
            << std::endl);
    }
    
    /*
     * Declare the pointers to the conservative variables and primitive variables.
     */
    
    std::vector<Real*> Q;
    Q.resize(d_num_eqn_conservative_var);
    
    std::vector<Real*> V;
    V.resize(d_num_eqn_primitive_var);
    
    int count_eqn = 0;
    
    /*
     * Convert primitive variables to conservative variables.
     */
    
    // Create the temporary side data.
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_density(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_pressure(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_internal_energy(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_primitive_var));
    
    Real* rho     = nullptr;
    Real* p       = nullptr;
    Real* epsilon = nullptr;
    
    /*
     * Get the thermodynamic properties of the species.
     */
    
    const int num_thermo_properties = d_equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<Real> thermo_properties;
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
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
        
        // Get the density.
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear index.
            const int idx_primitive_var = i + num_ghosts_0_primitive_var;
            
            rho[idx_primitive_var] = V[0][idx_primitive_var];
        }
        
        // Get the pressure.
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_primitive_var;
                i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                i++)
        {
            // Compute the linear index.
            const int idx_primitive_var = i + num_ghosts_0_primitive_var;
            
            p[idx_primitive_var] = V[1 + d_dim.getValue()][idx_primitive_var];
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->getEquationOfState()->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            thermo_properties_const_ptr,
            0);
        
        // Set the density.
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            
            Q[0][idx_conservative_var] = V[0][idx_primitive_var];
        }
        
        // Set the momentum.
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            
            Q[1][idx_conservative_var] = rho[idx_primitive_var]*V[1][idx_primitive_var];
        }
        
        // Set the total energy.
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_conservative_var;
                i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                i++)
        {
            // Compute the linear indices.
            const int idx_conservative_var = i + num_ghosts_0_conservative_var;
            const int idx_primitive_var    = i + num_ghosts_0_primitive_var;
            
            Q[2][idx_conservative_var] = rho[idx_primitive_var]*
                (epsilon[idx_primitive_var] + Real(1)/Real(2)*
                V[1][idx_primitive_var]*V[1][idx_primitive_var]);
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
        
        // Get the density.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_primitive_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                    i++)
            {
                // Compute the linear index.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                rho[idx_primitive_var] = V[0][idx_primitive_var];
            }
        }
        
        // Get the pressure.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_primitive_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                    i++)
            {
                // Compute the linear index.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                p[idx_primitive_var] = V[1 + d_dim.getValue()][idx_primitive_var];
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->getEquationOfState()->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            thermo_properties_const_ptr,
            0);
        
        // Set the density.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_conservative_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                Q[0][idx_conservative_var] = V[0][idx_primitive_var];
            }
        }
        
        // Set the momentum.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_conservative_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                Q[1][idx_conservative_var] = rho[idx_primitive_var]*V[1][idx_primitive_var];
                Q[2][idx_conservative_var] = rho[idx_primitive_var]*V[2][idx_primitive_var];
            }
        }
        
        // Set the total energy.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_conservative_var;
                    i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                    i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                Q[3][idx_conservative_var] = rho[idx_primitive_var]*
                    (epsilon[idx_primitive_var] + Real(1)/Real(2)*(
                    V[1][idx_primitive_var]*V[1][idx_primitive_var] +
                    V[2][idx_primitive_var]*V[2][idx_primitive_var]));
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
        
        // Get the density.
        for (int j = -num_ghosts_1_primitive_var;
                j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                rho[idx_primitive_var] = V[0][idx_primitive_var];
            }
        }
        
        // Get the pressure.
        for (int j = -num_ghosts_1_primitive_var;
                j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                p[idx_primitive_var] = V[1 + d_dim.getValue()][idx_primitive_var];
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->getEquationOfState()->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            thermo_properties_const_ptr,
            1);
        
        // Set the density.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                Q[0][idx_conservative_var] = V[0][idx_primitive_var];
            }
        }
        
        // Set the momentum.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                Q[1][idx_conservative_var] = rho[idx_primitive_var]*V[1][idx_primitive_var];
                Q[2][idx_conservative_var] = rho[idx_primitive_var]*V[2][idx_primitive_var];
            }
        }
        
        // Set the total energy.
        for (int j = -num_ghosts_1_primitive_var;
            j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
            j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_conservative_var = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                Q[3][idx_conservative_var] = rho[idx_primitive_var]*
                    (epsilon[idx_primitive_var] + Real(1)/Real(2)*(
                    V[1][idx_primitive_var]*V[1][idx_primitive_var] +
                    V[2][idx_primitive_var]*V[2][idx_primitive_var]));
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
        
        // Get the density.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_primitive_var;
                        i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                        i++)
                {
                    // Compute the linear index.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    rho[idx_primitive_var] = V[0][idx_primitive_var];
                }
            }
        }
        
        // Get the pressure.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_primitive_var;
                        i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                        i++)
                {
                    // Compute the linear index.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    p[idx_primitive_var] = V[1 + d_dim.getValue()][idx_primitive_var];
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->getEquationOfState()->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            thermo_properties_const_ptr,
            0);
        
        // Set the density.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    Q[0][idx_conservative_var] = V[0][idx_primitive_var];
                }
            }
        }
        
        // Set the momentum.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    Q[1][idx_conservative_var] = rho[idx_primitive_var]*V[1][idx_primitive_var];
                    Q[2][idx_conservative_var] = rho[idx_primitive_var]*V[2][idx_primitive_var];
                    Q[3][idx_conservative_var] = rho[idx_primitive_var]*V[3][idx_primitive_var];
                }
            }
        }
        
        // Set the total energy.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    Q[4][idx_conservative_var] = rho[idx_primitive_var]*
                        (epsilon[idx_primitive_var] + Real(1)/Real(2)*(
                        V[1][idx_primitive_var]*V[1][idx_primitive_var] + 
                        V[2][idx_primitive_var]*V[2][idx_primitive_var] +
                        V[3][idx_primitive_var]*V[3][idx_primitive_var]));
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
        
        // Get the density.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                    j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                    j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    rho[idx_primitive_var] = V[0][idx_primitive_var];
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
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    p[idx_primitive_var] = V[1 + d_dim.getValue()][idx_primitive_var];
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->getEquationOfState()->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            thermo_properties_const_ptr,
            1);
        
        // Set the density.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                    j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                    j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    Q[0][idx_conservative_var] = V[0][idx_primitive_var];
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
                HAMERS_PRAGMA_SIMD
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
                    
                    Q[1][idx_conservative_var] = rho[idx_primitive_var]*V[1][idx_primitive_var];
                    Q[2][idx_conservative_var] = rho[idx_primitive_var]*V[2][idx_primitive_var];
                    Q[3][idx_conservative_var] = rho[idx_primitive_var]*V[3][idx_primitive_var];
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
                HAMERS_PRAGMA_SIMD
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
                    
                    Q[4][idx_conservative_var] = rho[idx_primitive_var]*
                        (epsilon[idx_primitive_var] + Real(1)/Real(2)*(
                        V[1][idx_primitive_var]*V[1][idx_primitive_var] + 
                        V[2][idx_primitive_var]*V[2][idx_primitive_var] +
                        V[3][idx_primitive_var]*V[3][idx_primitive_var]));
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
        
        // Get the density.
        for (int si = 0; si < d_num_species; si++)
        {
            for (int k = -num_ghosts_2_primitive_var;
                 k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                            (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                            (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                                ghostcell_dim_1_primitive_var;
                        
                        rho[idx_primitive_var] = V[0][idx_primitive_var];
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
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_primitive_var = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    p[idx_primitive_var] = V[1 + d_dim.getValue()][idx_primitive_var];
                }
            }
        }
        
        // Compute the specific internal energy.
        d_equation_of_state_mixing_rules->getEquationOfState()->computeInternalEnergy(
            data_internal_energy,
            data_density,
            data_pressure,
            thermo_properties_const_ptr,
            2);
        
        // Set the density.
        for (int k = -num_ghosts_2_conservative_var;
                k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
                k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    Q[0][idx_conservative_var] = V[0][idx_primitive_var];
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
                HAMERS_PRAGMA_SIMD
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
                    
                    Q[1][idx_conservative_var] = rho[idx_primitive_var]*V[1][idx_primitive_var];
                    Q[2][idx_conservative_var] = rho[idx_primitive_var]*V[2][idx_primitive_var];
                    Q[3][idx_conservative_var] = rho[idx_primitive_var]*V[3][idx_primitive_var];
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
                HAMERS_PRAGMA_SIMD
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
                    
                    Q[4][idx_conservative_var] = rho[idx_primitive_var]*
                        (epsilon[idx_primitive_var] + Real(1)/Real(2)*(
                        V[1][idx_primitive_var]*V[1][idx_primitive_var] + 
                        V[2][idx_primitive_var]*V[2][idx_primitive_var] +
                        V[3][idx_primitive_var]*V[3][idx_primitive_var]));
                }
            }
        }
    }
}


/*
 * Check whether the given cell conservative variables are within the bounds.
 */
void
FlowModelBasicUtilitiesSingleSpecies::checkCellDataOfConservativeVariablesBounded(
    HAMERS_SHARED_PTR<pdat::CellData<int> >& bounded_flag,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_variables)
{
    NULL_USE(bounded_flag);
    NULL_USE(conservative_variables);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelBasicUtilitiesSingleSpecies::"
        << "checkCellDataOfConservativeVariablesBounded()\n"
        << "Method checkCellDataOfConservativeVariablesBounded()"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Check whether the given side conservative variables are within the bounds.
 */
void
FlowModelBasicUtilitiesSingleSpecies::checkSideDataOfConservativeVariablesBounded(
    HAMERS_SHARED_PTR<pdat::SideData<int> >& bounded_flag,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables)
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
            << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "checkSideDataOfConservativeVariablesBounded()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
                << std::endl);
        }
    }
    const hier::IntVector interior_dims_flag = bounded_flag->getBox().numberCells();
    if (interior_dims_flag != interior_dims)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "checkSideDataOfConservativeVariablesBounded()\n"
            << "The interior dimension of the flag does not match that of patch."
            << std::endl);
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_conservative_var != conservative_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "checkSideDataOfConservativeVariablesBounded()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_flag != num_ghosts_conservative_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "checkSideDataOfConservativeVariablesBounded()\n"
            << "The ghost cell width of the flag does not match that of conservative variables."
            << std::endl);
    }
    
    bounded_flag->fillAll(1);
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<Real*> Q;
    Q.resize(d_num_eqn);
    
    int* are_bounded = nullptr;
    
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
        
        // Check if density and total energy are bounded.
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_conservative_var;
             i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_conservative_var;
            
            if (Q[0][idx_face] > Real(0))
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
            
            if (Q[d_num_eqn - 1][idx_face] > Real(0))
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
        
        // Check if density and total energy are bounded.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_conservative_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1);
                
                if (Q[0][idx_face] > Real(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (Q[d_num_eqn - 1][idx_face] > Real(0))
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
        
        // Check if density and total energy are bounded.
        for (int j = -num_ghosts_1_conservative_var;
             j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
             j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_conservative_var) +
                    (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var;
                
                if (Q[0][idx_face] > Real(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (Q[d_num_eqn - 1][idx_face] > Real(0))
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
        
        // Check if density and total energy are bounded.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_conservative_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_conservative_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*(ghostcell_dim_0_conservative_var + 1) +
                        (k + num_ghosts_2_conservative_var)*(ghostcell_dim_0_conservative_var + 1)*
                            ghostcell_dim_1_conservative_var;
                    
                    if (Q[0][idx_face] > Real(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_eqn - 1][idx_face] > Real(0))
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
        
        // Check if density and total energy are bounded.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_conservative_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_conservative_var;
                 j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            (ghostcell_dim_1_conservative_var + 1);
                    
                    if (Q[0][idx_face] > Real(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_eqn - 1][idx_face] > Real(0))
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
        
        // Check if density and total energy are bounded.
        for (int k = -num_ghosts_2_conservative_var;
             k < interior_dim_2 + 1 + num_ghosts_2_conservative_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_conservative_var) +
                        (j + num_ghosts_1_conservative_var)*ghostcell_dim_0_conservative_var +
                        (k + num_ghosts_2_conservative_var)*ghostcell_dim_0_conservative_var*
                            ghostcell_dim_1_conservative_var;
                    
                    if (Q[0][idx_face] > Real(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (Q[d_num_eqn - 1][idx_face] > Real(0))
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
FlowModelBasicUtilitiesSingleSpecies::checkCellDataOfPrimitiveVariablesBounded(
    HAMERS_SHARED_PTR<pdat::CellData<int> >& bounded_flag,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& primitive_variables)
{
    NULL_USE(bounded_flag);
    NULL_USE(primitive_variables);
    
    TBOX_ERROR(d_object_name
        << ": FlowModelBasicUtilitiesSingleSpecies::"
        << "checkCellDataOfPrimitiveVariablesBounded()\n"
        << "Method checkCellDataOfPrimitiveVariablesBounded()"
        << " is not yet implemented."
        << std::endl);
}


/*
 * Check whether the given side primitive variables are within the bounds.
 */
void
FlowModelBasicUtilitiesSingleSpecies::checkSideDataOfPrimitiveVariablesBounded(
    HAMERS_SHARED_PTR<pdat::SideData<int> >& bounded_flag,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables)
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
            << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "checkSideDataOfPrimitiveVariablesBounded()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    const hier::IntVector interior_dims_flag = bounded_flag->getBox().numberCells();
    if (interior_dims_flag != interior_dims)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "checkSideDataOfPrimitiveVariablesBounded()\n"
            << "The interior dimension of the flag does not match that of patch."
            << std::endl);
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_primitive_var != primitive_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "checkSideDataOfPrimitiveVariablesBounded()\n"
                << "The primitive variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_flag != num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "checkSideDataOfPrimitiveVariablesBounded()\n"
            << "The ghost cell width of the flag does not match that of primitive variables."
            << std::endl);
    }
    
    bounded_flag->fillAll(1);
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<Real*> V;
    V.resize(d_num_eqn);
    
    int* are_bounded = nullptr;
    
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
        
        // Check if density and pressure are bounded.
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_primitive_var;
             i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_primitive_var;
            
            if (V[0][idx_face] > Real(0))
            {
                are_bounded[idx_face] &= 1;
            }
            else
            {
                are_bounded[idx_face] &= 0;
            }
            
            if (V[d_num_eqn - 1][idx_face] > Real(0))
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
        
        // Check if density and pressure are bounded.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_primitive_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1);
                
                if (V[0][idx_face] > Real(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (V[d_num_eqn - 1][idx_face] > Real(0))
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
        
        // Check if density and pressure are bounded.
        for (int j = -num_ghosts_1_primitive_var;
             j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
             j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_primitive_var) +
                    (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var;
                
                if (V[0][idx_face] > Real(0))
                {
                    are_bounded[idx_face] &= 1;
                }
                else
                {
                    are_bounded[idx_face] &= 0;
                }
                
                if (V[d_num_eqn - 1][idx_face] > Real(0))
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
        
        // Check if density and pressure are bounded.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_primitive_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_primitive_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*(ghostcell_dim_0_primitive_var + 1) +
                        (k + num_ghosts_2_primitive_var)*(ghostcell_dim_0_primitive_var + 1)*
                            ghostcell_dim_1_primitive_var;
                    
                    if (V[0][idx_face] > Real(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_eqn - 1][idx_face] > Real(0))
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
        
        // Check if density and pressure are bounded.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_primitive_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_primitive_var;
                 j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            (ghostcell_dim_1_primitive_var + 1);
                    
                    if (V[0][idx_face] > Real(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_eqn - 1][idx_face] > Real(0))
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
        
        // Check if density and pressure are bounded.
        for (int k = -num_ghosts_2_primitive_var;
             k < interior_dim_2 + 1 + num_ghosts_2_primitive_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_primitive_var) +
                        (j + num_ghosts_1_primitive_var)*ghostcell_dim_0_primitive_var +
                        (k + num_ghosts_2_primitive_var)*ghostcell_dim_0_primitive_var*
                            ghostcell_dim_1_primitive_var;
                    
                    if (V[0][idx_face] > Real(0))
                    {
                        are_bounded[idx_face] &= 1;
                    }
                    else
                    {
                        are_bounded[idx_face] &= 0;
                    }
                    
                    if (V[d_num_eqn - 1][idx_face] > Real(0))
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
FlowModelBasicUtilitiesSingleSpecies::registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables(
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
            << ": FlowModelBasicUtilitiesSingleSpecies::"
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
FlowModelBasicUtilitiesSingleSpecies::registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
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
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables()\n"
            << "No patch is registered yet."
            << std::endl);
    }
    
    d_proj_var_primitive_averaging_type = averaging_type;
    
    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
    
    num_subghosts_of_data.insert(
        std::pair<std::string, hier::IntVector>("SOUND_SPEED", num_subghosts));
    
    flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
}


/*
 * Get the number of projection variables for transformation between conservative
 * variables and characteristic variables.
 */
int
FlowModelBasicUtilitiesSingleSpecies::getNumberOfProjectionVariablesForConservativeVariables() const
{
    return d_dim.getValue() + 5;
}


/*
 * Get the number of projection variables for transformation between primitive variables
 * and characteristic variables.
 */
int
FlowModelBasicUtilitiesSingleSpecies::getNumberOfProjectionVariablesForPrimitiveVariables() const
{
    return 2;
}


/*
 * Compute the side data of the projection variables for transformation between conservative variables and
 * characteristic variables.
 */
void
FlowModelBasicUtilitiesSingleSpecies::computeSideDataOfProjectionVariablesForConservativeVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables)
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
    
    if (static_cast<int>(projection_variables.size()) != d_dim.getValue() + 5)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
            << "The number of projection variables are incorrect."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int vi = 0; vi < d_dim.getValue() + 5; vi++)
    {
        const hier::IntVector interior_dims_projection_var =
            projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int vi = 1; vi < d_dim.getValue() + 5; vi++)
    {
        if (num_ghosts_projection_var != projection_variables[vi]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                << "The projection variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var > num_ghosts)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
            << "The projection variables have ghost cell width larger than that of density."
            << std::endl);
    }
    
    // Get the cell data of the conservative variables.
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_density =
        flow_model_tmp->getCellData("DENSITY");
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_momentum =
        flow_model_tmp->getCellData("MOMENTUM");
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_total_energy =
        flow_model_tmp->getCellData("TOTAL_ENERGY");
    
    // Get the pointers to the cell data of density and total energy.
    
    Real* rho = data_density->getPointer(0);
    Real* E   = data_total_energy->getPointer(0);
    
    /*
     * Create temporary side data of averaged conservative variables.
     */
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_density_averaged(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_projection_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_momentum_averaged(
        new pdat::SideData<Real>(interior_box, d_dim.getValue(), num_ghosts_projection_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_total_energy_averaged(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_projection_var));
    
    /*
     * Create other temporary side data.
     */
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_internal_energy_averaged(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_projection_var));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > data_pressure_averaged(
        new pdat::SideData<Real>(interior_box, 1, num_ghosts_projection_var));
    
    /*
     * Declare pointers to the side data of averaged density and total energy.
     */
    
    Real* rho_average     = nullptr;
    Real* E_average       = nullptr;
    Real* e_average       = nullptr;
    Real* epsilon_average = nullptr;
    Real* p_average       = nullptr;
    Real* H_average       = nullptr;
    
    /*
     * Get the thermodynamic properties of the species.
     */
    
    const int num_thermo_properties = d_equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<Real> thermo_properties;
    std::vector<Real*> thermo_properties_ptr;
    std::vector<const Real*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_0_projection_var = num_ghosts_projection_var[0];
        
        // Get the pointer to the cell data of momentum.
        
        Real* rho_u = data_momentum->getPointer(0);
        
        // Declare pointers to the side data of averaged momentum and velocity.
        
        Real* rho_u_average = nullptr;
        Real* u_average     = nullptr;
        
        switch (d_proj_var_conservative_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the averaged conservative variables in the x-direction.
                 */
                
                rho_average   = data_density_averaged->getPointer(0);
                rho_u_average = data_momentum_averaged->getPointer(0, 0);
                E_average     = data_total_energy_averaged->getPointer(0);
                
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_projection_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i + num_ghosts_0_projection_var;
                    const int idx_L = i - 1 + num_ghosts_0;
                    const int idx_R = i + num_ghosts_0;
                    
                    rho_average[idx_face_x]   = Real(1)/Real(2)*(rho[idx_L] + rho[idx_R]);
                    rho_u_average[idx_face_x] = Real(1)/Real(2)*(rho_u[idx_L] + rho_u[idx_R]);
                    E_average[idx_face_x]     = Real(1)/Real(2)*(E[idx_L] + E[idx_R]);
                }
                
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                u_average       = projection_variables[0]->getPointer(0);
                e_average       = projection_variables[1]->getPointer(0);
                epsilon_average = data_internal_energy_averaged->getPointer(0);
                p_average       = data_pressure_averaged->getPointer(0);
                H_average       = projection_variables[2]->getPointer(0);
                
                // Compute the velocity and internal energy.
                
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_projection_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i + num_ghosts_0_projection_var;
                    
                    u_average[idx_face_x] = rho_u_average[idx_face_x]/rho_average[idx_face_x];
                    
                    e_average[idx_face_x] = E_average[idx_face_x]/rho_average[idx_face_x];
                    
                    epsilon_average[idx_face_x] = e_average[idx_face_x] -
                        Real(1)/Real(2)*u_average[idx_face_x]*u_average[idx_face_x];
                }
                
                // Compute the presure, sound speed, partial pressure partial density and Gruneisen parameter.
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                    data_pressure_averaged,
                    data_density_averaged,
                    data_internal_energy_averaged,
                    thermo_properties_const_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                    projection_variables[3],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressureDerivativeWithDensity(
                    projection_variables[4],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                    projection_variables[5],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    0);
                
                // Compute the total specific enthalpy.
                
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_projection_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i + num_ghosts_0_projection_var;
                    
                    H_average[idx_face_x] = e_average[idx_face_x] +
                        p_average[idx_face_x]/rho_average[idx_face_x];
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                    << "Unknown d_proj_var_conservative_averaging_type given."
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
        
        // Get the pointers to the cell data of momentum.
        
        Real* rho_u = data_momentum->getPointer(0);
        Real* rho_v = data_momentum->getPointer(1);
        
        // Declare pointers to the side data of averaged momentum and velocity.
        
        Real* rho_u_average = nullptr;
        Real* rho_v_average = nullptr;
        Real* u_average     = nullptr;
        Real* v_average     = nullptr;
        
        switch (d_proj_var_conservative_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the averaged conservative variables in the x-direction.
                 */
                
                rho_average   = data_density_averaged->getPointer(0);
                rho_u_average = data_momentum_averaged->getPointer(0, 0);
                rho_v_average = data_momentum_averaged->getPointer(0, 1);
                E_average     = data_total_energy_averaged->getPointer(0);
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
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
                        
                        rho_average[idx_face_x]   = Real(1)/Real(2)*(rho[idx_L] + rho[idx_R]);
                        rho_u_average[idx_face_x] = Real(1)/Real(2)*(rho_u[idx_L] + rho_u[idx_R]);
                        rho_v_average[idx_face_x] = Real(1)/Real(2)*(rho_v[idx_L] + rho_v[idx_R]);
                        E_average[idx_face_x]     = Real(1)/Real(2)*(E[idx_L] + E[idx_R]);
                    }
                }
                
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                u_average       = projection_variables[0]->getPointer(0);
                v_average       = projection_variables[1]->getPointer(0);
                e_average       = projection_variables[2]->getPointer(0);
                epsilon_average = data_internal_energy_averaged->getPointer(0);
                p_average       = data_pressure_averaged->getPointer(0);
                H_average       = projection_variables[3]->getPointer(0);
                
                // Compute the velocity and internal energy.
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = -num_ghosts_0_projection_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = (i + num_ghosts_0_projection_var) +
                            (j + num_ghosts_1_projection_var)*(ghostcell_dim_0_projection_var + 1);
                        
                        u_average[idx_face_x] = rho_u_average[idx_face_x]/rho_average[idx_face_x];
                        v_average[idx_face_x] = rho_v_average[idx_face_x]/rho_average[idx_face_x];
                        
                        e_average[idx_face_x] = E_average[idx_face_x]/rho_average[idx_face_x];
                        
                        epsilon_average[idx_face_x] = e_average[idx_face_x] -
                            Real(1)/Real(2)*(u_average[idx_face_x]*u_average[idx_face_x] +
                                v_average[idx_face_x]*v_average[idx_face_x]);
                    }
                }
                
                // Compute the presure, sound speed, partial pressure partial density and Gruneisen parameter.
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                    data_pressure_averaged,
                    data_density_averaged,
                    data_internal_energy_averaged,
                    thermo_properties_const_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                    projection_variables[4],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressureDerivativeWithDensity(
                    projection_variables[5],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                    projection_variables[6],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    0);
                
                // Compute the total specific enthalpy.
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = -num_ghosts_0_projection_var;
                         i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                         i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = (i + num_ghosts_0_projection_var) +
                            (j + num_ghosts_1_projection_var)*(ghostcell_dim_0_projection_var + 1);
                        
                        H_average[idx_face_x] = e_average[idx_face_x] +
                            p_average[idx_face_x]/rho_average[idx_face_x];
                    }
                }
                
                /*
                 * Compute the averaged conservative variables in the y-direction.
                 */
                
                rho_average   = data_density_averaged->getPointer(1);
                rho_u_average = data_momentum_averaged->getPointer(1, 0);
                rho_v_average = data_momentum_averaged->getPointer(1, 1);
                E_average     = data_total_energy_averaged->getPointer(1);
                
                for (int j = -num_ghosts_1_projection_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                     j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = (i + num_ghosts_0_projection_var) +
                            (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var;
                        
                        const int idx_B = (i + num_ghosts_0) +
                            (j - 1 + num_ghosts_1)*ghostcell_dim_0;
                        
                        const int idx_T = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0;
                        
                        rho_average[idx_face_y]   = Real(1)/Real(2)*(rho[idx_B] + rho[idx_T]);
                        rho_u_average[idx_face_y] = Real(1)/Real(2)*(rho_u[idx_B] + rho_u[idx_T]);
                        rho_v_average[idx_face_y] = Real(1)/Real(2)*(rho_v[idx_B] + rho_v[idx_T]);
                        E_average[idx_face_y]     = Real(1)/Real(2)*(E[idx_B] + E[idx_T]);
                    }
                }
                
                /*
                 * Compute the projection variables in the y-direction.
                 */
                
                u_average       = projection_variables[0]->getPointer(1);
                v_average       = projection_variables[1]->getPointer(1);
                e_average       = projection_variables[2]->getPointer(1);
                epsilon_average = data_internal_energy_averaged->getPointer(1);
                p_average       = data_pressure_averaged->getPointer(1);
                H_average       = projection_variables[3]->getPointer(1);
                
                // Compute the velocity and internal energy.
                
                for (int j = -num_ghosts_1_projection_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                     j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = (i + num_ghosts_0_projection_var) +
                            (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var;
                        
                        u_average[idx_face_y] = rho_u_average[idx_face_y]/rho_average[idx_face_y];
                        v_average[idx_face_y] = rho_v_average[idx_face_y]/rho_average[idx_face_y];
                        
                        e_average[idx_face_y] = E_average[idx_face_y]/rho_average[idx_face_y];
                        
                        epsilon_average[idx_face_y] = e_average[idx_face_y] -
                            Real(1)/Real(2)*(u_average[idx_face_y]*u_average[idx_face_y] +
                                v_average[idx_face_y]*v_average[idx_face_y]);
                    }
                }
                
                // Compute the presure, sound speed, partial pressure partial density and Gruneisen parameter.
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                    data_pressure_averaged,
                    data_density_averaged,
                    data_internal_energy_averaged,
                    thermo_properties_const_ptr,
                    1);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                    projection_variables[4],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    1);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressureDerivativeWithDensity(
                    projection_variables[5],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    1);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                    projection_variables[6],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    1);
                
                // Compute the total specific enthalpy.
                
                for (int j = -num_ghosts_1_projection_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                     j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = (i + num_ghosts_0_projection_var) +
                            (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var;
                        
                        H_average[idx_face_y] = e_average[idx_face_y] +
                            p_average[idx_face_y]/rho_average[idx_face_y];
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                    << "Unknown d_proj_var_conservative_averaging_type given."
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
        
        // Get the pointers to the cell data of momentum.
        
        Real* rho_u = data_momentum->getPointer(0);
        Real* rho_v = data_momentum->getPointer(1);
        Real* rho_w = data_momentum->getPointer(2);
        
        // Declare pointers to the side data of averaged momentum and velocity.
        
        Real* rho_u_average = nullptr;
        Real* rho_v_average = nullptr;
        Real* rho_w_average = nullptr;
        Real* u_average     = nullptr;
        Real* v_average     = nullptr;
        Real* w_average     = nullptr;
        
        switch (d_proj_var_conservative_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the averaged conservative variables in the x-direction.
                 */
                
                rho_average   = data_density_averaged->getPointer(0);
                rho_u_average = data_momentum_averaged->getPointer(0, 0);
                rho_v_average = data_momentum_averaged->getPointer(0, 1);
                rho_w_average = data_momentum_averaged->getPointer(0, 2);
                E_average     = data_total_energy_averaged->getPointer(0);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                            
                            rho_average[idx_face_x]   = Real(1)/Real(2)*(rho[idx_L] + rho[idx_R]);
                            rho_u_average[idx_face_x] = Real(1)/Real(2)*(rho_u[idx_L] + rho_u[idx_R]);
                            rho_v_average[idx_face_x] = Real(1)/Real(2)*(rho_v[idx_L] + rho_v[idx_R]);
                            rho_w_average[idx_face_x] = Real(1)/Real(2)*(rho_w[idx_L] + rho_w[idx_R]);
                            E_average[idx_face_x]     = Real(1)/Real(2)*(E[idx_L] + E[idx_R]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                u_average       = projection_variables[0]->getPointer(0);
                v_average       = projection_variables[1]->getPointer(0);
                w_average       = projection_variables[2]->getPointer(0);
                e_average       = projection_variables[3]->getPointer(0);
                epsilon_average = data_internal_energy_averaged->getPointer(0);
                p_average       = data_pressure_averaged->getPointer(0);
                H_average       = projection_variables[4]->getPointer(0);
                
                // Compute the velocity and internal energy.
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = -num_ghosts_0_projection_var;
                             i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*(ghostcell_dim_0_projection_var + 1) +
                                (k + num_ghosts_2_projection_var)*(ghostcell_dim_0_projection_var + 1)*
                                    ghostcell_dim_1_projection_var;
                            
                            u_average[idx_face_x] = rho_u_average[idx_face_x]/rho_average[idx_face_x];
                            v_average[idx_face_x] = rho_v_average[idx_face_x]/rho_average[idx_face_x];
                            w_average[idx_face_x] = rho_w_average[idx_face_x]/rho_average[idx_face_x];
                            
                            e_average[idx_face_x] = E_average[idx_face_x]/rho_average[idx_face_x];
                            
                            epsilon_average[idx_face_x] = e_average[idx_face_x] -
                                Real(1)/Real(2)*(u_average[idx_face_x]*u_average[idx_face_x] +
                                    v_average[idx_face_x]*v_average[idx_face_x] +
                                    w_average[idx_face_x]*w_average[idx_face_x]);
                        }
                    }
                }
                
                // Compute the presure, sound speed, partial pressure partial density and Gruneisen parameter.
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                    data_pressure_averaged,
                    data_density_averaged,
                    data_internal_energy_averaged,
                    thermo_properties_const_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                    projection_variables[5],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressureDerivativeWithDensity(
                    projection_variables[6],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    0);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                    projection_variables[7],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    0);
                
                // Compute the total specific enthalpy.
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = -num_ghosts_0_projection_var;
                             i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                             i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*(ghostcell_dim_0_projection_var + 1) +
                                (k + num_ghosts_2_projection_var)*(ghostcell_dim_0_projection_var + 1)*
                                    ghostcell_dim_1_projection_var;
                            
                            H_average[idx_face_x] = e_average[idx_face_x] +
                                p_average[idx_face_x]/rho_average[idx_face_x];
                        }
                    }
                }
                
                /*
                 * Compute the averaged conservative variables in the y-direction.
                 */
                
                rho_average   = data_density_averaged->getPointer(1);
                rho_u_average = data_momentum_averaged->getPointer(1, 0);
                rho_v_average = data_momentum_averaged->getPointer(1, 1);
                rho_w_average = data_momentum_averaged->getPointer(1, 2);
                E_average     = data_total_energy_averaged->getPointer(1);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = -num_ghosts_1_projection_var;
                         j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                         j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                            
                            rho_average[idx_face_y]   = Real(1)/Real(2)*(rho[idx_B] + rho[idx_T]);
                            rho_u_average[idx_face_y] = Real(1)/Real(2)*(rho_u[idx_B] + rho_u[idx_T]);
                            rho_v_average[idx_face_y] = Real(1)/Real(2)*(rho_v[idx_B] + rho_v[idx_T]);
                            rho_w_average[idx_face_y] = Real(1)/Real(2)*(rho_w[idx_B] + rho_w[idx_T]);
                            E_average[idx_face_y]     = Real(1)/Real(2)*(E[idx_B] + E[idx_T]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the y-direction.
                 */
                
                u_average       = projection_variables[0]->getPointer(1);
                v_average       = projection_variables[1]->getPointer(1);
                w_average       = projection_variables[2]->getPointer(1);
                e_average       = projection_variables[3]->getPointer(1);
                epsilon_average = data_internal_energy_averaged->getPointer(1);
                p_average       = data_pressure_averaged->getPointer(1);
                H_average       = projection_variables[4]->getPointer(1);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = -num_ghosts_1_projection_var;
                         j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                         j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var +
                                (k + num_ghosts_2_projection_var)*ghostcell_dim_0_projection_var*
                                    (ghostcell_dim_1_projection_var + 1);
                            
                            u_average[idx_face_y] = rho_u_average[idx_face_y]/rho_average[idx_face_y];
                            v_average[idx_face_y] = rho_v_average[idx_face_y]/rho_average[idx_face_y];
                            w_average[idx_face_y] = rho_w_average[idx_face_y]/rho_average[idx_face_y];
                            
                            e_average[idx_face_y] = E_average[idx_face_y]/rho_average[idx_face_y];
                            
                            epsilon_average[idx_face_y] = e_average[idx_face_y] -
                                Real(1)/Real(2)*(u_average[idx_face_y]*u_average[idx_face_y] +
                                    v_average[idx_face_y]*v_average[idx_face_y] +
                                    w_average[idx_face_y]*w_average[idx_face_y]);
                        }
                    }
                }
                
                // Compute the presure, sound speed, partial pressure partial density and Gruneisen parameter.
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                    data_pressure_averaged,
                    data_density_averaged,
                    data_internal_energy_averaged,
                    thermo_properties_const_ptr,
                    1);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                    projection_variables[5],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    1);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressureDerivativeWithDensity(
                    projection_variables[6],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    1);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                    projection_variables[7],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    1);
                
                // Compute the total specific enthalpy.
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = -num_ghosts_1_projection_var;
                         j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                         j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var +
                                (k + num_ghosts_2_projection_var)*ghostcell_dim_0_projection_var*
                                    (ghostcell_dim_1_projection_var + 1);
                            
                            H_average[idx_face_y] = e_average[idx_face_y] +
                                p_average[idx_face_y]/rho_average[idx_face_y];
                        }
                    }
                }
                
                /*
                 * Compute the averaged conservative variables in the z-direction.
                 */
                
                rho_average   = data_density_averaged->getPointer(2);
                rho_u_average = data_momentum_averaged->getPointer(2, 0);
                rho_v_average = data_momentum_averaged->getPointer(2, 1);
                rho_w_average = data_momentum_averaged->getPointer(2, 2);
                E_average     = data_total_energy_averaged->getPointer(2);
                
                for (int k = -num_ghosts_2_projection_var;
                     k < interior_dim_2 + 1 + num_ghosts_2_projection_var;
                     k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                            
                            rho_average[idx_face_z]   = Real(1)/Real(2)*(rho[idx_B] + rho[idx_F]);
                            rho_u_average[idx_face_z] = Real(1)/Real(2)*(rho_u[idx_B] + rho_u[idx_F]);
                            rho_v_average[idx_face_z] = Real(1)/Real(2)*(rho_v[idx_B] + rho_v[idx_F]);
                            rho_w_average[idx_face_z] = Real(1)/Real(2)*(rho_w[idx_B] + rho_w[idx_F]);
                            E_average[idx_face_z]     = Real(1)/Real(2)*(E[idx_B] + E[idx_F]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the z-direction.
                 */
                
                u_average       = projection_variables[0]->getPointer(2);
                v_average       = projection_variables[1]->getPointer(2);
                w_average       = projection_variables[2]->getPointer(2);
                e_average       = projection_variables[3]->getPointer(2);
                epsilon_average = data_internal_energy_averaged->getPointer(2);
                p_average       = data_pressure_averaged->getPointer(2);
                H_average       = projection_variables[4]->getPointer(2);
                
                for (int k = -num_ghosts_2_projection_var;
                     k < interior_dim_2 + 1 + num_ghosts_2_projection_var;
                     k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var +
                                (k + num_ghosts_2_projection_var)*ghostcell_dim_0_projection_var*
                                    ghostcell_dim_1_projection_var;
                            
                            u_average[idx_face_z] = rho_u_average[idx_face_z]/rho_average[idx_face_z];
                            v_average[idx_face_z] = rho_v_average[idx_face_z]/rho_average[idx_face_z];
                            w_average[idx_face_z] = rho_w_average[idx_face_z]/rho_average[idx_face_z];
                            
                            e_average[idx_face_z] = E_average[idx_face_z]/rho_average[idx_face_z];
                            
                            epsilon_average[idx_face_z] = e_average[idx_face_z] -
                                Real(1)/Real(2)*(u_average[idx_face_z]*u_average[idx_face_z] +
                                    v_average[idx_face_z]*v_average[idx_face_z] +
                                    w_average[idx_face_z]*w_average[idx_face_z]);
                        }
                    }
                }
                
                // Compute the presure, sound speed, partial pressure partial density and Gruneisen parameter.
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressure(
                    data_pressure_averaged,
                    data_density_averaged,
                    data_internal_energy_averaged,
                    thermo_properties_const_ptr,
                    2);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeSoundSpeed(
                    projection_variables[5],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    2);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computePressureDerivativeWithDensity(
                    projection_variables[6],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    2);
                
                d_equation_of_state_mixing_rules->getEquationOfState()->computeGruneisenParameter(
                    projection_variables[7],
                    data_density_averaged,
                    data_pressure_averaged,
                    thermo_properties_const_ptr,
                    2);
                
                // Compute the total specific enthalpy.
                
                for (int k = -num_ghosts_2_projection_var;
                     k < interior_dim_2 + 1 + num_ghosts_2_projection_var;
                     k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = (i + num_ghosts_0_projection_var) +
                                (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var +
                                (k + num_ghosts_2_projection_var)*ghostcell_dim_0_projection_var*
                                    ghostcell_dim_1_projection_var;
                            
                            H_average[idx_face_z] = e_average[idx_face_z] +
                                p_average[idx_face_z]/rho_average[idx_face_z];
                        }
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForConservativeVariables()\n"
                    << "Unknown d_proj_var_conservative_averaging_type given."
                    << std::endl);
            }
        }
    }
}


/*
 * Compute the side data of the projection variables for transformation between primitive variables and
 * characteristic variables.
 */
void
FlowModelBasicUtilitiesSingleSpecies::computeSideDataOfProjectionVariablesForPrimitiveVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables)
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
    
    if (static_cast<int>(projection_variables.size()) != 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
            << "There should be two projection variables."
            << std::endl);
    }
    
    /*
     * Check potential failures.
     */
    
    for (int vi = 0; vi < 2; vi++)
    {
        const hier::IntVector interior_dims_projection_var =
            projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var != projection_variables[1]->getGhostCellWidth())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables don't have same ghost cell width."
            << std::endl);
    }
    
    if (num_ghosts_projection_var > num_ghosts)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables have ghost cell width larger than that of density."
            << std::endl);
    }
    
    // Get the cell data of the variable density.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_density =
        flow_model_tmp->getCellData("DENSITY");
    
    // Get the cell data of sound speed.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > data_sound_speed =
        flow_model_tmp->getCellData("SOUND_SPEED");
    
    // Get the number of ghost cells and ghost cell dimension of sound speed.
    const hier::IntVector& num_subghosts_sound_speed = data_sound_speed->getGhostCellWidth();
    const hier::IntVector subghostcell_dims_sound_speed = data_sound_speed->getGhostBox().numberCells();
    
    if (num_ghosts_projection_var > num_subghosts_sound_speed)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
            << "The projection variables have ghost cell width larger than that of sound speed."
            << std::endl);
    }
    
    // Get the pointers to the cell data of density and sound speed.
    Real* rho = data_density->getPointer(0);
    Real* c = data_sound_speed->getPointer(0);
    
    /*
     * Declare pointers to different data.
     */
    
    Real* rho_average = nullptr;
    Real* c_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_0_projection_var = num_ghosts_projection_var[0];
        const int num_subghosts_0_sound_speed = num_subghosts_sound_speed[0];
        
        switch (d_proj_var_primitive_averaging_type)
        {
            case AVERAGING::SIMPLE:
            {
                /*
                 * Compute the projection variables in the x-direction.
                 */
                
                rho_average = projection_variables[0]->getPointer(0);
                c_average = projection_variables[1]->getPointer(0);
                
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_projection_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_projection_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i + num_ghosts_0_projection_var;
                    const int idx_L = i - 1 + num_ghosts_0;
                    const int idx_R = i + num_ghosts_0;
                    const int idx_sound_speed_L = i - 1 + num_subghosts_0_sound_speed;
                    const int idx_sound_speed_R = i + num_subghosts_0_sound_speed;
                    
                    rho_average[idx_face_x] = Real(1)/Real(2)*(rho[idx_L] + rho[idx_R]);
                    c_average[idx_face_x] = Real(1)/Real(2)*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                
                rho_average = projection_variables[0]->getPointer(0);
                c_average = projection_variables[1]->getPointer(0);
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
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
                        
                        const int idx_sound_speed_L = (i - 1 + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_sound_speed_R = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        rho_average[idx_face_x] = Real(1)/Real(2)*(rho[idx_L] + rho[idx_R]);
                        c_average[idx_face_x] = Real(1)/Real(2)*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                    }
                }
                
                /*
                 * Compute the projection variables in the y-direction.
                 */
                rho_average = projection_variables[0]->getPointer(1);
                c_average = projection_variables[1]->getPointer(1);
                
                for (int j = -num_ghosts_1_projection_var;
                     j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                     j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = (i + num_ghosts_0_projection_var) +
                            (j + num_ghosts_1_projection_var)*ghostcell_dim_0_projection_var;
                        
                        const int idx_B = (i + num_ghosts_0) +
                            (j - 1 + num_ghosts_1)*ghostcell_dim_0;
                        
                        const int idx_T = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0;
                        
                        const int idx_sound_speed_B = (i + num_subghosts_0_sound_speed) +
                            (j - 1 + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_sound_speed_T = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        rho_average[idx_face_y] = Real(1)/Real(2)*(rho[idx_B] + rho[idx_T]);
                        c_average[idx_face_y] = Real(1)/Real(2)*(c[idx_sound_speed_B] + c[idx_sound_speed_T]);
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                
                rho_average = projection_variables[0]->getPointer(0);
                c_average = projection_variables[1]->getPointer(0);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                            
                            const int idx_sound_speed_L = (i - 1 + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_sound_speed_R = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            rho_average[idx_face_x] = Real(1)/Real(2)*(rho[idx_L] + rho[idx_R]);
                            c_average[idx_face_x] = Real(1)/Real(2)*(c[idx_sound_speed_L] + c[idx_sound_speed_R]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the y-direction.
                 */
                
                rho_average = projection_variables[0]->getPointer(1);
                c_average = projection_variables[1]->getPointer(1);
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = -num_ghosts_1_projection_var;
                         j < interior_dim_1 + 1 + num_ghosts_1_projection_var;
                         j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                            
                            const int idx_sound_speed_B = (i + num_subghosts_0_sound_speed) +
                                (j - 1 + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_sound_speed_T = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            rho_average[idx_face_y] = Real(1)/Real(2)*(rho[idx_B] + rho[idx_T]);
                            c_average[idx_face_y] = Real(1)/Real(2)*(c[idx_sound_speed_B] + c[idx_sound_speed_T]);
                        }
                    }
                }
                
                /*
                 * Compute the projection variables in the z-direction.
                 */
                
                rho_average = projection_variables[0]->getPointer(2);
                c_average = projection_variables[1]->getPointer(2);
                
                for (int k = -num_ghosts_2_projection_var;
                     k < interior_dim_2 + 1 + num_ghosts_2_projection_var;
                     k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                            
                            const int idx_sound_speed_B = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k - 1 + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_sound_speed_F = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            rho_average[idx_face_z] = Real(1)/Real(2)*(rho[idx_B] + rho[idx_F]);
                            c_average[idx_face_z] = Real(1)/Real(2)*(c[idx_sound_speed_B] + c[idx_sound_speed_F]);
                        }
                    }
                }
                
                break;
            }
            case AVERAGING::ROE:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesSingleSpecies::"
                    << "computeSideDataOfProjectionVariablesForPrimitiveVariables()\n"
                    << "Roe averaging is not yet implemented."
                    << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBasicUtilitiesSingleSpecies::"
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
FlowModelBasicUtilitiesSingleSpecies::computeSideDataOfCharacteristicVariablesFromConservativeVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& characteristic_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables,
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
    
    std::vector<hier::IntVector> num_ghosts_conservative_var;
    num_ghosts_conservative_var.reserve(static_cast<int>(conservative_variables.size()));
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        num_ghosts_conservative_var.push_back(conservative_variables[vi]->
            getGhostCellWidth());
    }
    
    const hier::IntVector num_ghosts_projection_var = projection_variables[0]->getGhostCellWidth();
    
    /*
     * Get the ghost cell dimensions of characteristic and conservative variables.
     */
    
    const hier::IntVector ghostcell_dims_characteristic_var = characteristic_variables[0]->
        getGhostBox().numberCells();
    
    std::vector<hier::IntVector> ghostcell_dims_conservative_var;
    ghostcell_dims_conservative_var.reserve(static_cast<int>(conservative_variables.size()));
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        ghostcell_dims_conservative_var.push_back(conservative_variables[vi]->
            getGhostBox().numberCells());
    }
    
    /*
     * Check the size of variables.
     */
    
    if (static_cast<int>(characteristic_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(conservative_variables.size()) != 3)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    if (conservative_variables[0]->getDepth() != 1 ||
        conservative_variables[1]->getDepth() != d_dim.getValue() ||
        conservative_variables[2]->getDepth() != 1)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
            << "The depths of one or more conservative variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != d_dim.getValue() + 5)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
            << "The number of projection variables are incorrect."
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
                << "The interior dimension of the characteristic variables does not match that of patch."
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < d_dim.getValue() + 5; vi++)
    {
        const hier::IntVector interior_dims_projection_var = projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_characteristic_var != characteristic_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < d_dim.getValue() + 5; ei++)
    {
        if (num_ghosts_projection_var != projection_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
                << "The projection variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
            << "The ghost cell width of the projection variables does not match that of"
            << " characteristic variables."
            << std::endl);
    }
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        if (num_ghosts_conservative_var[vi] - num_ghosts_characteristic_var
                + (hier::IntVector::getOne(d_dim))*idx_offset < hier::IntVector::getZero(d_dim) ||
            num_ghosts_characteristic_var - num_ghosts_conservative_var[vi]
                + (hier::IntVector::getOne(d_dim))*(idx_offset + 1) > hier::IntVector::getZero(d_dim))
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromConservativeVariables()\n"
                << "The offset index is too large or the number of ghost of characteristic variable"
                << " is too large."
                << std::endl);
        }
    }
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<Real*> Q;
    Q.reserve(d_num_eqn);
    
    for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
    {
        int depth = conservative_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            Q.push_back(conservative_variables[vi]->getPointer(di));
        }
    }
    
    std::vector<Real*> W;
    W.resize(d_num_eqn);
    
    Real* e_average     = nullptr;
    Real* c_average     = nullptr;
    Real* Psi_average   = nullptr;
    Real* Gamma_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_0_rho = num_ghosts_conservative_var[0][0];
        const int num_ghosts_0_mom = num_ghosts_conservative_var[1][0];
        const int num_ghosts_0_E = num_ghosts_conservative_var[2][0];
        
        const int idx_offset_rho = idx_offset;
        const int idx_offset_mom = idx_offset;
        const int idx_offset_E = idx_offset;
        
        // Declare pointer to side data of velocity.
        
        Real* u_average = nullptr;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        u_average     = projection_variables[0]->getPointer(0);
        e_average     = projection_variables[1]->getPointer(0);
        c_average     = projection_variables[3]->getPointer(0);
        Psi_average   = projection_variables[4]->getPointer(0);
        Gamma_average = projection_variables[5]->getPointer(0);
        
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_characteristic_var;
             i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
             i++)
        {
            // Compute the linear indices.
            const int idx_face = i + num_ghosts_0_characteristic_var;
            const int idx_rho = i + idx_offset_rho + num_ghosts_0_rho;
            const int idx_mom = i + idx_offset_mom + num_ghosts_0_mom;
            const int idx_E = i + idx_offset_E + num_ghosts_0_E;
            
            W[0][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] -
                e_average[idx_face]) + Psi_average[idx_face] + u_average[idx_face]*c_average[idx_face])*
                Q[0][idx_rho] -
                (Gamma_average[idx_face]*u_average[idx_face] + c_average[idx_face])*Q[1][idx_mom] +
                Gamma_average[idx_face]*Q[2][idx_E])/
                    (Real(2)*c_average[idx_face]*c_average[idx_face]);
            
            W[1][idx_face] = (-(Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] -
                e_average[idx_face]) + Psi_average[idx_face] - c_average[idx_face]*c_average[idx_face])*
                Q[0][idx_rho] +
                Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] - Gamma_average[idx_face]*Q[2][idx_E])/
                    (c_average[idx_face]*c_average[idx_face]);
            
            W[2][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] -
                e_average[idx_face]) + Psi_average[idx_face] - u_average[idx_face]*c_average[idx_face])*
                Q[0][idx_rho] -
                (Gamma_average[idx_face]*u_average[idx_face] - c_average[idx_face])*Q[1][idx_mom] +
                Gamma_average[idx_face]*Q[2][idx_E])/
                    (Real(2)*c_average[idx_face]*c_average[idx_face]);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        
        const int num_ghosts_0_rho = num_ghosts_conservative_var[0][0];
        const int num_ghosts_1_rho = num_ghosts_conservative_var[0][1];
        const int ghostcell_dim_0_rho = ghostcell_dims_conservative_var[0][0];
        
        const int num_ghosts_0_mom = num_ghosts_conservative_var[1][0];
        const int num_ghosts_1_mom = num_ghosts_conservative_var[1][1];
        const int ghostcell_dim_0_mom = ghostcell_dims_conservative_var[1][0];
        
        const int num_ghosts_0_E = num_ghosts_conservative_var[2][0];
        const int num_ghosts_1_E = num_ghosts_conservative_var[2][1];
        const int ghostcell_dim_0_E = ghostcell_dims_conservative_var[2][0];
        
        const int idx_offset_rho = idx_offset;
        const int idx_offset_mom = idx_offset;
        const int idx_offset_E = idx_offset;
        
        // Declare pointers to side data of velocity.
        
        Real* u_average = nullptr;
        Real* v_average = nullptr;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        u_average     = projection_variables[0]->getPointer(0);
        v_average     = projection_variables[1]->getPointer(0);
        e_average     = projection_variables[2]->getPointer(0);
        c_average     = projection_variables[4]->getPointer(0);
        Psi_average   = projection_variables[5]->getPointer(0);
        Gamma_average = projection_variables[6]->getPointer(0);
        
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_characteristic_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                 i++)
            {
                // Compute the linear indices.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1);
                
                const int idx_rho = (i + idx_offset_rho + num_ghosts_0_rho) +
                    (j + num_ghosts_1_rho)*ghostcell_dim_0_rho;
                
                const int idx_mom = (i + idx_offset_mom + num_ghosts_0_mom) +
                    (j + num_ghosts_1_mom)*ghostcell_dim_0_mom;
                
                const int idx_E = (i + idx_offset_E + num_ghosts_0_E) +
                    (j + num_ghosts_1_E)*ghostcell_dim_0_E;
                
                W[0][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                    v_average[idx_face]*v_average[idx_face] - e_average[idx_face]) + Psi_average[idx_face] +
                    u_average[idx_face]*c_average[idx_face])*Q[0][idx_rho] -
                    (Gamma_average[idx_face]*u_average[idx_face] + c_average[idx_face])*Q[1][idx_mom] -
                    Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] +
                    Gamma_average[idx_face]*Q[3][idx_E])/
                        (Real(2)*c_average[idx_face]*c_average[idx_face]);
                
                W[1][idx_face] = (-(Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                    v_average[idx_face]*v_average[idx_face] - e_average[idx_face]) + Psi_average[idx_face] -
                    c_average[idx_face]*c_average[idx_face])*Q[0][idx_rho] +
                    Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] +
                    Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] -
                    Gamma_average[idx_face]*Q[3][idx_E])/
                        (c_average[idx_face]*c_average[idx_face]);
                
                W[2][idx_face] = -Q[0][idx_rho] + Q[2][idx_mom]/v_average[idx_face];
                
                W[3][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                    v_average[idx_face]*v_average[idx_face] - e_average[idx_face]) + Psi_average[idx_face] -
                    u_average[idx_face]*c_average[idx_face])*Q[0][idx_rho] -
                    (Gamma_average[idx_face]*u_average[idx_face] - c_average[idx_face])*Q[1][idx_mom] -
                    Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] +
                    Gamma_average[idx_face]*Q[3][idx_E])/
                        (Real(2)*c_average[idx_face]*c_average[idx_face]);
            }
        }
        
        /*
         * Compute the characteristic variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        u_average     = projection_variables[0]->getPointer(1);
        v_average     = projection_variables[1]->getPointer(1);
        e_average     = projection_variables[2]->getPointer(1);
        c_average     = projection_variables[4]->getPointer(1);
        Psi_average   = projection_variables[5]->getPointer(1);
        Gamma_average = projection_variables[6]->getPointer(1);
        
        for (int j = -num_ghosts_1_characteristic_var;
             j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
             j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var;
                
                const int idx_rho = (i + num_ghosts_0_rho) +
                    (j + idx_offset_rho + num_ghosts_1_rho)*ghostcell_dim_0_rho;
                
                const int idx_mom = (i + num_ghosts_0_mom) +
                    (j + idx_offset_mom + num_ghosts_1_mom)*ghostcell_dim_0_mom;
                
                const int idx_E = (i + num_ghosts_0_E) +
                    (j + idx_offset_E + num_ghosts_1_E)*ghostcell_dim_0_E;
                
                W[0][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                    v_average[idx_face]*v_average[idx_face] - e_average[idx_face]) + Psi_average[idx_face] +
                    v_average[idx_face]*c_average[idx_face])*Q[0][idx_rho] -
                    Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] -
                    (Gamma_average[idx_face]*v_average[idx_face] + c_average[idx_face])*Q[2][idx_mom] +
                    Gamma_average[idx_face]*Q[3][idx_E])/
                        (Real(2)*c_average[idx_face]*c_average[idx_face]);
                
                W[1][idx_face] = -Q[0][idx_rho] + Q[1][idx_mom]/u_average[idx_face];
                
                W[2][idx_face] = (-(Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                    v_average[idx_face]*v_average[idx_face] - e_average[idx_face]) + Psi_average[idx_face] -
                    c_average[idx_face]*c_average[idx_face])*Q[0][idx_rho] +
                    Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] +
                    Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] -
                    Gamma_average[idx_face]*Q[3][idx_E])/
                        (c_average[idx_face]*c_average[idx_face]);
                
                W[3][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                    v_average[idx_face]*v_average[idx_face] - e_average[idx_face]) + Psi_average[idx_face] -
                    v_average[idx_face]*c_average[idx_face])*Q[0][idx_rho] -
                    Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] -
                    (Gamma_average[idx_face]*v_average[idx_face] - c_average[idx_face])*Q[2][idx_mom] +
                    Gamma_average[idx_face]*Q[3][idx_E])/
                        (Real(2)*c_average[idx_face]*c_average[idx_face]);
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
        
        const int num_ghosts_0_rho = num_ghosts_conservative_var[0][0];
        const int num_ghosts_1_rho = num_ghosts_conservative_var[0][1];
        const int num_ghosts_2_rho = num_ghosts_conservative_var[0][2];
        const int ghostcell_dim_0_rho = ghostcell_dims_conservative_var[0][0];
        const int ghostcell_dim_1_rho = ghostcell_dims_conservative_var[0][1];
        
        const int num_ghosts_0_mom = num_ghosts_conservative_var[1][0];
        const int num_ghosts_1_mom = num_ghosts_conservative_var[1][1];
        const int num_ghosts_2_mom = num_ghosts_conservative_var[1][2];
        const int ghostcell_dim_0_mom = ghostcell_dims_conservative_var[1][0];
        const int ghostcell_dim_1_mom = ghostcell_dims_conservative_var[1][1];
        
        const int num_ghosts_0_E = num_ghosts_conservative_var[2][0];
        const int num_ghosts_1_E = num_ghosts_conservative_var[2][1];
        const int num_ghosts_2_E = num_ghosts_conservative_var[2][2];
        const int ghostcell_dim_0_E = ghostcell_dims_conservative_var[2][0];
        const int ghostcell_dim_1_E = ghostcell_dims_conservative_var[2][1];
        
        const int idx_offset_rho = idx_offset;
        const int idx_offset_mom = idx_offset;
        const int idx_offset_E = idx_offset;
        
        // Declare pointers to side data of velocity.
        
        Real* u_average = nullptr;
        Real* v_average = nullptr;
        Real* w_average = nullptr;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        u_average     = projection_variables[0]->getPointer(0);
        v_average     = projection_variables[1]->getPointer(0);
        w_average     = projection_variables[2]->getPointer(0);
        e_average     = projection_variables[3]->getPointer(0);
        c_average     = projection_variables[5]->getPointer(0);
        Psi_average   = projection_variables[6]->getPointer(0);
        Gamma_average = projection_variables[7]->getPointer(0);
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_characteristic_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1) +
                        (k + num_ghosts_2_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1)*
                            ghostcell_dim_1_characteristic_var;
                    
                    const int idx_rho = (i + idx_offset_rho + num_ghosts_0_rho) +
                        (j + num_ghosts_1_rho)*ghostcell_dim_0_rho +
                        (k + num_ghosts_2_rho)*ghostcell_dim_0_rho*
                            ghostcell_dim_1_rho;
                    
                    const int idx_mom = (i + idx_offset_mom + num_ghosts_0_mom) +
                        (j + num_ghosts_1_mom)*ghostcell_dim_0_mom +
                        (k + num_ghosts_2_mom)*ghostcell_dim_0_mom*
                            ghostcell_dim_1_mom;
                    
                    const int idx_E = (i + idx_offset_E + num_ghosts_0_E) +
                        (j + num_ghosts_1_E)*ghostcell_dim_0_E +
                        (k + num_ghosts_2_E)*ghostcell_dim_0_E*
                            ghostcell_dim_1_E;
                    
                    W[0][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] + u_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] -
                        (Gamma_average[idx_face]*u_average[idx_face] + c_average[idx_face])*Q[1][idx_mom] -
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] -
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] +
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (Real(2)*c_average[idx_face]*c_average[idx_face]);
                    
                    W[1][idx_face] = (-(Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] - c_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] +
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] +
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] +
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] -
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (c_average[idx_face]*c_average[idx_face]);
                    
                    W[2][idx_face] = -Q[0][idx_rho] + Q[2][idx_mom]/v_average[idx_face];
                    
                    W[3][idx_face] = -Q[0][idx_rho] + Q[3][idx_mom]/w_average[idx_face];
                    
                    W[4][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] - u_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] -
                        (Gamma_average[idx_face]*u_average[idx_face] - c_average[idx_face])*Q[1][idx_mom] -
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] -
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] +
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (Real(2)*c_average[idx_face]*c_average[idx_face]);
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
        
        u_average     = projection_variables[0]->getPointer(1);
        v_average     = projection_variables[1]->getPointer(1);
        w_average     = projection_variables[2]->getPointer(1);
        e_average     = projection_variables[3]->getPointer(1);
        c_average     = projection_variables[5]->getPointer(1);
        Psi_average   = projection_variables[6]->getPointer(1);
        Gamma_average = projection_variables[7]->getPointer(1);
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_characteristic_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                 j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            (ghostcell_dim_1_characteristic_var + 1);
                    
                    const int idx_rho = (i + num_ghosts_0_rho) +
                        (j + idx_offset_rho + num_ghosts_1_rho)*ghostcell_dim_0_rho +
                        (k + num_ghosts_2_rho)*ghostcell_dim_0_rho*
                            ghostcell_dim_1_rho;
                    
                    const int idx_mom = (i + num_ghosts_0_mom) +
                        (j + idx_offset_mom + num_ghosts_1_mom)*ghostcell_dim_0_mom +
                        (k + num_ghosts_2_mom)*ghostcell_dim_0_mom*
                            ghostcell_dim_1_mom;
                    
                    const int idx_E = (i + num_ghosts_0_E) +
                        (j + idx_offset_E + num_ghosts_1_E)*ghostcell_dim_0_E +
                        (k + num_ghosts_2_E)*ghostcell_dim_0_E*
                            ghostcell_dim_1_E;
                    
                    W[0][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] + w_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] -
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] -
                        (Gamma_average[idx_face]*v_average[idx_face] + c_average[idx_face])*Q[2][idx_mom] -
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] +
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (Real(2)*c_average[idx_face]*c_average[idx_face]);
                    
                    W[1][idx_face] = -Q[0][idx_rho] + Q[1][idx_mom]/u_average[idx_face];
                    
                    W[2][idx_face] = (-(Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] - c_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] +
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] +
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] +
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] -
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (c_average[idx_face]*c_average[idx_face]);
                    
                    W[3][idx_face] = -Q[0][idx_rho] + Q[3][idx_mom]/w_average[idx_face];
                    
                    W[4][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] - w_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] -
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] -
                        (Gamma_average[idx_face]*v_average[idx_face] - c_average[idx_face])*Q[2][idx_mom] -
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] +
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (Real(2)*c_average[idx_face]*c_average[idx_face]);
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
        
        u_average     = projection_variables[0]->getPointer(2);
        v_average     = projection_variables[1]->getPointer(2);
        w_average     = projection_variables[2]->getPointer(2);
        e_average     = projection_variables[3]->getPointer(2);
        c_average     = projection_variables[5]->getPointer(2);
        Psi_average   = projection_variables[6]->getPointer(2);
        Gamma_average = projection_variables[7]->getPointer(2);
        
        for (int k = -num_ghosts_2_characteristic_var;
             k < interior_dim_2 + 1 + num_ghosts_2_characteristic_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            ghostcell_dim_1_characteristic_var;
                    
                    const int idx_rho = (i + num_ghosts_0_rho) +
                        (j + num_ghosts_1_rho)*ghostcell_dim_0_rho +
                        (k + idx_offset_rho + num_ghosts_2_rho)*ghostcell_dim_0_rho*
                            ghostcell_dim_1_rho;
                    
                    const int idx_mom = (i + num_ghosts_0_mom) +
                        (j + num_ghosts_1_mom)*ghostcell_dim_0_mom +
                        (k + idx_offset_mom + num_ghosts_2_mom)*ghostcell_dim_0_mom*
                            ghostcell_dim_1_mom;
                    
                    const int idx_E = (i + num_ghosts_0_E) +
                        (j + num_ghosts_1_E)*ghostcell_dim_0_E +
                        (k + idx_offset_E + num_ghosts_2_E)*ghostcell_dim_0_E*
                            ghostcell_dim_1_E;
                    
                    W[0][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] + w_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] -
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] -
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] -
                        (Gamma_average[idx_face]*w_average[idx_face] + c_average[idx_face])*Q[3][idx_mom] +
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (Real(2)*c_average[idx_face]*c_average[idx_face]);
                    
                    W[1][idx_face] = -Q[0][idx_rho] + Q[1][idx_mom]/u_average[idx_face];
                    
                    W[2][idx_face] = -Q[0][idx_rho] + Q[2][idx_mom]/v_average[idx_face];
                    
                    W[3][idx_face] = (-(Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] - c_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] +
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] +
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] +
                        Gamma_average[idx_face]*w_average[idx_face]*Q[3][idx_mom] -
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (c_average[idx_face]*c_average[idx_face]);
                    
                    W[4][idx_face] = ((Gamma_average[idx_face]*(u_average[idx_face]*u_average[idx_face] +
                        v_average[idx_face]*v_average[idx_face] + w_average[idx_face]*w_average[idx_face] -
                        e_average[idx_face]) + Psi_average[idx_face] - w_average[idx_face]*c_average[idx_face])*
                        Q[0][idx_rho] -
                        Gamma_average[idx_face]*u_average[idx_face]*Q[1][idx_mom] -
                        Gamma_average[idx_face]*v_average[idx_face]*Q[2][idx_mom] -
                        (Gamma_average[idx_face]*w_average[idx_face] - c_average[idx_face])*Q[3][idx_mom] +
                        Gamma_average[idx_face]*Q[4][idx_E])/
                            (Real(2)*c_average[idx_face]*c_average[idx_face]);
                }
            }
        }
    }
}


/*
 * Compute the side data of characteristic variables from primitive variables.
 */
void
FlowModelBasicUtilitiesSingleSpecies::computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& characteristic_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& primitive_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables,
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
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(primitive_variables.size()) != 3)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    if (primitive_variables[0]->getDepth() != 1 ||
        primitive_variables[1]->getDepth() != d_dim.getValue() ||
        primitive_variables[2]->getDepth() != 1)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The depths of one or more primitive variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "There should be two projection variables."
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The interior dimension of the primitive variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < 2; vi++)
    {
        const hier::IntVector interior_dims_projection_var = projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    if (num_ghosts_projection_var != projection_variables[1]->getGhostCellWidth())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
            << "The projection variables don't have same ghost cell width."
            << std::endl);
    }
    
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfCharacteristicVariablesFromPrimitiveVariables()\n"
                << "The offset index is too large or the number of ghost of characteristic variable"
                << " is too large."
                << std::endl);
        }
    }
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<Real*> V;
    V.reserve(d_num_eqn);
    
    for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
    {
        int depth = primitive_variables[vi]->getDepth();
        
        for (int di = 0; di < depth; di++)
        {
            V.push_back(primitive_variables[vi]->getPointer(di));
        }
    }
    
    std::vector<Real*> W;
    W.resize(d_num_eqn);
    
    Real* rho_average = nullptr;
    Real* c_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_0_rho = num_ghosts_primitive_var[0][0];
        const int num_ghosts_0_vel = num_ghosts_primitive_var[1][0];
        const int num_ghosts_0_p = num_ghosts_primitive_var[2][0];
        
        const int idx_offset_rho = idx_offset;
        const int idx_offset_vel = idx_offset;
        const int idx_offset_p = idx_offset;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        rho_average = projection_variables[0]->getPointer(0);
        c_average = projection_variables[1]->getPointer(0);
        
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_characteristic_var;
             i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
             i++)
        {
            // Compute the linear indices.
            const int idx_face = i + num_ghosts_0_characteristic_var;
            const int idx_rho = i + idx_offset_rho + num_ghosts_0_rho;
            const int idx_vel = i + idx_offset_vel + num_ghosts_0_vel;
            const int idx_p = i + idx_offset_p + num_ghosts_0_p;
            
            W[0][idx_face] = -Real(1)/Real(2)*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] +
                Real(1)/Real(2)*V[2][idx_p];
            W[1][idx_face] = V[0][idx_rho] - Real(1)/(c_average[idx_face]*c_average[idx_face])*V[2][idx_p];
            W[2][idx_face] = Real(1)/Real(2)*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] +
                Real(1)/Real(2)*V[2][idx_p];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        
        const int num_ghosts_0_rho = num_ghosts_primitive_var[0][0];
        const int num_ghosts_1_rho = num_ghosts_primitive_var[0][1];
        const int ghostcell_dim_0_rho = ghostcell_dims_primitive_var[0][0];
        
        const int num_ghosts_0_vel = num_ghosts_primitive_var[1][0];
        const int num_ghosts_1_vel = num_ghosts_primitive_var[1][1];
        const int ghostcell_dim_0_vel = ghostcell_dims_primitive_var[1][0];
        
        const int num_ghosts_0_p = num_ghosts_primitive_var[2][0];
        const int num_ghosts_1_p = num_ghosts_primitive_var[2][1];
        const int ghostcell_dim_0_p = ghostcell_dims_primitive_var[2][0];
        
        const int idx_offset_rho = idx_offset;
        const int idx_offset_vel = idx_offset;
        const int idx_offset_p = idx_offset;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        rho_average = projection_variables[0]->getPointer(0);
        c_average = projection_variables[1]->getPointer(0);
        
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_characteristic_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                 i++)
            {
                // Compute the linear indices.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1);
                
                const int idx_rho = (i + idx_offset_rho + num_ghosts_0_rho) +
                    (j + num_ghosts_1_rho)*ghostcell_dim_0_rho;
                
                const int idx_vel = (i + idx_offset_vel + num_ghosts_0_vel) +
                    (j + num_ghosts_1_vel)*ghostcell_dim_0_vel;
                
                const int idx_p = (i + idx_offset_p + num_ghosts_0_p) +
                    (j + num_ghosts_1_p)*ghostcell_dim_0_p;
                
                W[0][idx_face] = -Real(1)/Real(2)*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] +
                    Real(1)/Real(2)*V[3][idx_p];
                W[1][idx_face] = V[0][idx_rho] - Real(1)/(c_average[idx_face]*c_average[idx_face])*V[3][idx_p];
                W[2][idx_face] = V[2][idx_vel];
                W[3][idx_face] = Real(1)/Real(2)*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] +
                    Real(1)/Real(2)*V[3][idx_p];
            }
        }
        
        /*
         * Compute the characteristic variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        rho_average = projection_variables[0]->getPointer(1);
        c_average = projection_variables[1]->getPointer(1);
        
        for (int j = -num_ghosts_1_characteristic_var;
             j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
             j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var;
                
                const int idx_rho = (i + num_ghosts_0_rho) +
                    (j + idx_offset_rho + num_ghosts_1_rho)*ghostcell_dim_0_rho;
                
                const int idx_vel = (i + num_ghosts_0_vel) +
                    (j + idx_offset_vel + num_ghosts_1_vel)*ghostcell_dim_0_vel;
                
                const int idx_p = (i + num_ghosts_0_p) +
                    (j + idx_offset_p + num_ghosts_1_p)*ghostcell_dim_0_p;
                
                W[0][idx_face] = -Real(1)/Real(2)*rho_average[idx_face]*c_average[idx_face]*V[2][idx_vel] +
                    Real(1)/Real(2)*V[3][idx_p];
                W[1][idx_face] = V[0][idx_rho] - Real(1)/(c_average[idx_face]*c_average[idx_face])*V[3][idx_p];
                W[2][idx_face] = V[1][idx_vel];
                W[3][idx_face] = Real(1)/Real(2)*rho_average[idx_face]*c_average[idx_face]*V[2][idx_vel] +
                    Real(1)/Real(2)*V[3][idx_p];
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
        
        const int num_ghosts_0_rho = num_ghosts_primitive_var[0][0];
        const int num_ghosts_1_rho = num_ghosts_primitive_var[0][1];
        const int num_ghosts_2_rho = num_ghosts_primitive_var[0][2];
        const int ghostcell_dim_0_rho = ghostcell_dims_primitive_var[0][0];
        const int ghostcell_dim_1_rho = ghostcell_dims_primitive_var[0][1];
        
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
        
        const int idx_offset_rho = idx_offset;
        const int idx_offset_vel = idx_offset;
        const int idx_offset_p = idx_offset;
        
        /*
         * Compute the characteristic variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        rho_average = projection_variables[0]->getPointer(0);
        c_average = projection_variables[1]->getPointer(0);
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_characteristic_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                     i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1) +
                        (k + num_ghosts_2_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1)*
                            ghostcell_dim_1_characteristic_var;
                    
                    const int idx_rho = (i + idx_offset_rho + num_ghosts_0_rho) +
                        (j + num_ghosts_1_rho)*ghostcell_dim_0_rho +
                        (k + num_ghosts_2_rho)*ghostcell_dim_0_rho*
                            ghostcell_dim_1_rho;
                    
                    const int idx_vel = (i + idx_offset_vel + num_ghosts_0_vel) +
                        (j + num_ghosts_1_vel)*ghostcell_dim_0_vel +
                        (k + num_ghosts_2_vel)*ghostcell_dim_0_vel*
                            ghostcell_dim_1_vel;
                    
                    const int idx_p = (i + idx_offset_p + num_ghosts_0_p) +
                        (j + num_ghosts_1_p)*ghostcell_dim_0_p +
                        (k + num_ghosts_2_p)*ghostcell_dim_0_p*
                            ghostcell_dim_1_p;
                    
                    W[0][idx_face] = -Real(1)/Real(2)*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] +
                        Real(1)/Real(2)*V[4][idx_p];
                    W[1][idx_face] = V[0][idx_rho] - Real(1)/(c_average[idx_face]*c_average[idx_face])*V[4][idx_p];
                    W[2][idx_face] = V[2][idx_vel];
                    W[3][idx_face] = V[3][idx_vel];
                    W[4][idx_face] = Real(1)/Real(2)*rho_average[idx_face]*c_average[idx_face]*V[1][idx_vel] +
                        Real(1)/Real(2)*V[4][idx_p];
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
        
        rho_average = projection_variables[0]->getPointer(1);
        c_average = projection_variables[1]->getPointer(1);
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_characteristic_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                 j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            (ghostcell_dim_1_characteristic_var + 1);
                    
                    const int idx_rho = (i + num_ghosts_0_rho) +
                        (j + idx_offset_rho + num_ghosts_1_rho)*ghostcell_dim_0_rho +
                        (k + num_ghosts_2_rho)*ghostcell_dim_0_rho*
                            ghostcell_dim_1_rho;
                    
                    const int idx_vel = (i + num_ghosts_0_vel) +
                        (j + idx_offset_vel + num_ghosts_1_vel)*ghostcell_dim_0_vel +
                        (k + num_ghosts_2_vel)*ghostcell_dim_0_vel*
                            ghostcell_dim_1_vel;
                    
                    const int idx_p = (i + num_ghosts_0_p) +
                        (j + idx_offset_p + num_ghosts_1_p)*ghostcell_dim_0_p +
                        (k + num_ghosts_2_p)*ghostcell_dim_0_p*
                            ghostcell_dim_1_p;
                    
                    W[0][idx_face] = -Real(1)/Real(2)*rho_average[idx_face]*c_average[idx_face]*V[2][idx_vel] +
                        Real(1)/Real(2)*V[4][idx_p];
                    W[1][idx_face] = V[0][idx_rho] - Real(1)/(c_average[idx_face]*c_average[idx_face])*V[4][idx_p];
                    W[2][idx_face] = V[1][idx_vel];
                    W[3][idx_face] = V[3][idx_vel];
                    W[4][idx_face] = Real(1)/Real(2)*rho_average[idx_face]*c_average[idx_face]*V[2][idx_vel] +
                        Real(1)/Real(2)*V[4][idx_p];
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
        
        rho_average = projection_variables[0]->getPointer(2);
        c_average = projection_variables[1]->getPointer(2);
        
        for (int k = -num_ghosts_2_characteristic_var;
             k < interior_dim_2 + 1 + num_ghosts_2_characteristic_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            ghostcell_dim_1_characteristic_var;
                    
                    const int idx_rho = (i + num_ghosts_0_rho) +
                        (j + num_ghosts_1_rho)*ghostcell_dim_0_rho +
                        (k + idx_offset_rho + num_ghosts_2_rho)*ghostcell_dim_0_rho*
                            ghostcell_dim_1_rho;
                    
                    const int idx_vel = (i + num_ghosts_0_vel) +
                        (j + num_ghosts_1_vel)*ghostcell_dim_0_vel +
                        (k + idx_offset_vel + num_ghosts_2_vel)*ghostcell_dim_0_vel*
                            ghostcell_dim_1_vel;
                    
                    const int idx_p = (i + num_ghosts_0_p) +
                        (j + num_ghosts_1_p)*ghostcell_dim_0_p +
                        (k + idx_offset_p + num_ghosts_2_p)*ghostcell_dim_0_p*
                            ghostcell_dim_1_p;
                    
                    W[0][idx_face] = -Real(1)/Real(2)*rho_average[idx_face]*c_average[idx_face]*V[3][idx_vel] +
                        Real(1)/Real(2)*V[4][idx_p];
                    W[1][idx_face] = V[0][idx_rho] - Real(1)/(c_average[idx_face]*c_average[idx_face])*V[4][idx_p];
                    W[2][idx_face] = V[1][idx_vel];
                    W[3][idx_face] = V[2][idx_vel];
                    W[4][idx_face] = Real(1)/Real(2)*rho_average[idx_face]*c_average[idx_face]*V[3][idx_vel] +
                        Real(1)/Real(2)*V[4][idx_p];
                }
            }
        }
    }
}


/*
 * Compute the side data of conservative variables from characteristic variables.
 */
void
FlowModelBasicUtilitiesSingleSpecies::computeSideDataOfConservativeVariablesFromCharacteristicVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& characteristic_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables)
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
    
    const hier::IntVector num_ghosts_conservative_var = conservative_variables[0]->
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
    
    if (static_cast<int>(conservative_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(characteristic_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != d_dim.getValue() + 5)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
            << "The number of projection variables are incorrect."
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the conservative variables does not match that of patch."
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the characteristic variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < d_dim.getValue() + 5; vi++)
    {
        const hier::IntVector interior_dims_projection_var = projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the projection variables does not match that of patch."
                << std::endl);
        }
    }
    
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_conservative_var != conservative_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
                << "The conservative variables don't have same ghost cell width."
                << std::endl);
        }
    }
    for (int ei = 1; ei < d_num_eqn; ei++)
    {
        if (num_ghosts_characteristic_var != characteristic_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    for (int ei = 1; ei < d_dim.getValue() + 5; ei++)
    {
        if (num_ghosts_projection_var != projection_variables[ei]->getGhostCellWidth())
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
                << "The projection variables don't have same ghost cell width."
                << std::endl);
        }
    }
    
    if (num_ghosts_projection_var != num_ghosts_conservative_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of conservative"
            << " variables."
            << std::endl);
    }
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfConservativeVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of characteristic"
            << " variables."
            << std::endl);
    }
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<Real*> Q;
    std::vector<Real*> W;
    Q.resize(d_num_eqn);
    W.resize(d_num_eqn);
    
    Real* e_average     = nullptr;
    Real* H_average     = nullptr;
    Real* c_average     = nullptr;
    Real* Psi_average   = nullptr;
    Real* Gamma_average = nullptr;
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        
        // Declare pointer to side data of velocity.
        
        Real* u_average = nullptr;
        
        /*
         * Compute the conservative variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        u_average     = projection_variables[0]->getPointer(0);
        e_average     = projection_variables[1]->getPointer(0);
        H_average     = projection_variables[2]->getPointer(0);
        c_average     = projection_variables[3]->getPointer(0);
        Psi_average   = projection_variables[4]->getPointer(0);
        Gamma_average = projection_variables[5]->getPointer(0);
        
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_characteristic_var;
             i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_characteristic_var;
            
            Q[0][idx_face] = W[0][idx_face] + W[1][idx_face] + W[2][idx_face];
            
            Q[1][idx_face] = (u_average[idx_face] - c_average[idx_face])*W[0][idx_face] +
                u_average[idx_face]*W[1][idx_face] +
                (u_average[idx_face] + c_average[idx_face])*W[2][idx_face];
            
            Q[2][idx_face] = (H_average[idx_face] - u_average[idx_face]*c_average[idx_face])*
                W[0][idx_face] +
                (e_average[idx_face] - Psi_average[idx_face]/Gamma_average[idx_face])*W[1][idx_face] +
                (H_average[idx_face] + u_average[idx_face]*c_average[idx_face])*
                W[2][idx_face];
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_characteristic_var = num_ghosts_characteristic_var[0];
        const int num_ghosts_1_characteristic_var = num_ghosts_characteristic_var[1];
        const int ghostcell_dim_0_characteristic_var = ghostcell_dims_characteristic_var[0];
        
        // Declare pointers to side data of velocity.
        
        Real* u_average = nullptr;
        Real* v_average = nullptr;
        
        /*
         * Compute the conservative variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        u_average     = projection_variables[0]->getPointer(0);
        v_average     = projection_variables[1]->getPointer(0);
        e_average     = projection_variables[2]->getPointer(0);
        H_average     = projection_variables[3]->getPointer(0);
        c_average     = projection_variables[4]->getPointer(0);
        Psi_average   = projection_variables[5]->getPointer(0);
        Gamma_average = projection_variables[6]->getPointer(0);
        
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_characteristic_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1);
                
                Q[0][idx_face] = W[0][idx_face] + W[1][idx_face] + W[3][idx_face];
                
                Q[1][idx_face] = (u_average[idx_face] - c_average[idx_face])*W[0][idx_face] +
                    u_average[idx_face]*W[1][idx_face] +
                    (u_average[idx_face] + c_average[idx_face])*W[3][idx_face];
                
                Q[2][idx_face] = v_average[idx_face]*(W[0][idx_face] + W[1][idx_face] + W[2][idx_face] +
                    W[3][idx_face]);
                
                Q[3][idx_face] = (H_average[idx_face] - u_average[idx_face]*c_average[idx_face])*
                    W[0][idx_face] +
                    (e_average[idx_face] - Psi_average[idx_face]/Gamma_average[idx_face])*W[1][idx_face] +
                    v_average[idx_face]*v_average[idx_face]*W[2][idx_face] +
                    (H_average[idx_face] + u_average[idx_face]*c_average[idx_face])*
                    W[3][idx_face];
            }
        }
        
        /*
         * Compute the conservative variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(1);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        u_average     = projection_variables[0]->getPointer(1);
        v_average     = projection_variables[1]->getPointer(1);
        e_average     = projection_variables[2]->getPointer(1);
        H_average     = projection_variables[3]->getPointer(1);
        c_average     = projection_variables[4]->getPointer(1);
        Psi_average   = projection_variables[5]->getPointer(1);
        Gamma_average = projection_variables[6]->getPointer(1);
        
        for (int j = -num_ghosts_1_characteristic_var;
             j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
             j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var;
                
                Q[0][idx_face] = W[0][idx_face] + W[2][idx_face] + W[3][idx_face];
                
                Q[1][idx_face] = u_average[idx_face]*(W[0][idx_face] + W[1][idx_face] + W[2][idx_face] +
                    W[3][idx_face]);
                
                Q[2][idx_face] = (v_average[idx_face] - c_average[idx_face])*W[0][idx_face] +
                    v_average[idx_face]*W[2][idx_face] +
                    (v_average[idx_face] + c_average[idx_face])*W[3][idx_face];
                
                Q[3][idx_face] = (H_average[idx_face] - v_average[idx_face]*c_average[idx_face])*
                    W[0][idx_face] +
                    u_average[idx_face]*u_average[idx_face]*W[1][idx_face] +
                    (e_average[idx_face] - Psi_average[idx_face]/Gamma_average[idx_face])*W[2][idx_face] +
                    (H_average[idx_face] + v_average[idx_face]*c_average[idx_face])*
                    W[3][idx_face];
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
        
        // Declare pointers to side data of velocity.
        
        Real* u_average = nullptr;
        Real* v_average = nullptr;
        Real* w_average = nullptr;
        
        /*
         * Compute the conservative variables in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(0);
        }
        
        u_average     = projection_variables[0]->getPointer(0);
        v_average     = projection_variables[1]->getPointer(0);
        w_average     = projection_variables[2]->getPointer(0);
        e_average     = projection_variables[3]->getPointer(0);
        H_average     = projection_variables[4]->getPointer(0);
        c_average     = projection_variables[5]->getPointer(0);
        Psi_average   = projection_variables[6]->getPointer(0);
        Gamma_average = projection_variables[7]->getPointer(0);
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_characteristic_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1) +
                        (k + num_ghosts_2_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1)*
                            ghostcell_dim_1_characteristic_var;
                    
                    Q[0][idx_face] = W[0][idx_face] + W[1][idx_face] + W[4][idx_face];
                    
                    Q[1][idx_face] = (u_average[idx_face] - c_average[idx_face])*W[0][idx_face] +
                        u_average[idx_face]*W[1][idx_face] +
                        (u_average[idx_face] + c_average[idx_face])*W[4][idx_face];
                    
                    Q[2][idx_face] = v_average[idx_face]*(W[0][idx_face] + W[1][idx_face] + W[2][idx_face] +
                        W[4][idx_face]);
                    
                    Q[3][idx_face] = w_average[idx_face]*(W[0][idx_face] + W[1][idx_face] + W[3][idx_face] +
                        W[4][idx_face]);
                    
                    Q[4][idx_face] = (H_average[idx_face] - u_average[idx_face]*c_average[idx_face])*
                        W[0][idx_face] +
                        (e_average[idx_face] - Psi_average[idx_face]/Gamma_average[idx_face])*W[1][idx_face] +
                        v_average[idx_face]*v_average[idx_face]*W[2][idx_face] +
                        w_average[idx_face]*w_average[idx_face]*W[3][idx_face] +
                        (H_average[idx_face] + u_average[idx_face]*c_average[idx_face])*
                        W[4][idx_face];
                }
            }
        }
        
        /*
         * Compute the conservative variables in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(1);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(1);
        }
        
        u_average     = projection_variables[0]->getPointer(1);
        v_average     = projection_variables[1]->getPointer(1);
        w_average     = projection_variables[2]->getPointer(1);
        e_average     = projection_variables[3]->getPointer(1);
        H_average     = projection_variables[4]->getPointer(1);
        c_average     = projection_variables[5]->getPointer(1);
        Psi_average   = projection_variables[6]->getPointer(1);
        Gamma_average = projection_variables[7]->getPointer(1);
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_characteristic_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                 j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            (ghostcell_dim_1_characteristic_var + 1);
                    
                    Q[0][idx_face] = W[0][idx_face] + W[2][idx_face] + W[4][idx_face];
                    
                    Q[1][idx_face] = u_average[idx_face]*(W[0][idx_face] + W[1][idx_face] + W[2][idx_face] +
                        W[4][idx_face]);
                    
                    Q[2][idx_face] = (v_average[idx_face] - c_average[idx_face])*W[0][idx_face] +
                        v_average[idx_face]*W[2][idx_face] +
                        (v_average[idx_face] + c_average[idx_face])*W[4][idx_face];
                    
                    Q[3][idx_face] = w_average[idx_face]*(W[0][idx_face] + W[2][idx_face] + W[3][idx_face] +
                        W[4][idx_face]);
                    
                    Q[4][idx_face] = (H_average[idx_face] - v_average[idx_face]*c_average[idx_face])*
                        W[0][idx_face] +
                        u_average[idx_face]*u_average[idx_face]*W[1][idx_face] +
                        (e_average[idx_face] - Psi_average[idx_face]/Gamma_average[idx_face])*W[2][idx_face] +
                        w_average[idx_face]*w_average[idx_face]*W[3][idx_face] +
                        (H_average[idx_face] + v_average[idx_face]*c_average[idx_face])*
                        W[4][idx_face];
                }
            }
        }
        
        /*
         * Compute the conservative variables in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q[ei] = conservative_variables[ei]->getPointer(2);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            W[ei] = characteristic_variables[ei]->getPointer(2);
        }
        
        u_average     = projection_variables[0]->getPointer(2);
        v_average     = projection_variables[1]->getPointer(2);
        w_average     = projection_variables[2]->getPointer(2);
        e_average     = projection_variables[3]->getPointer(2);
        H_average     = projection_variables[4]->getPointer(2);
        c_average     = projection_variables[5]->getPointer(2);
        Psi_average   = projection_variables[6]->getPointer(2);
        Gamma_average = projection_variables[7]->getPointer(2);
        
        for (int k = -num_ghosts_2_characteristic_var;
             k < interior_dim_2 + 1 + num_ghosts_2_characteristic_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            ghostcell_dim_1_characteristic_var;
                    
                    Q[0][idx_face] = W[0][idx_face] + W[3][idx_face] + W[4][idx_face];
                    
                    Q[1][idx_face] = u_average[idx_face]*(W[0][idx_face] + W[1][idx_face] + W[3][idx_face] +
                        W[4][idx_face]);
                    
                    Q[2][idx_face] = v_average[idx_face]*(W[0][idx_face] + W[2][idx_face] + W[3][idx_face] +
                        W[4][idx_face]);
                    
                    Q[3][idx_face] = (w_average[idx_face] - c_average[idx_face])*W[0][idx_face] +
                        w_average[idx_face]*W[3][idx_face] +
                        (w_average[idx_face] + c_average[idx_face])*W[4][idx_face];
                    
                    Q[4][idx_face] = (H_average[idx_face] - w_average[idx_face]*c_average[idx_face])*
                        W[0][idx_face] +
                        u_average[idx_face]*u_average[idx_face]*W[1][idx_face] +
                        v_average[idx_face]*v_average[idx_face]*W[2][idx_face] +
                        (e_average[idx_face] - Psi_average[idx_face]/Gamma_average[idx_face])*W[3][idx_face] +
                        (H_average[idx_face] + w_average[idx_face]*c_average[idx_face])*
                        W[4][idx_face];
                }
            }
        }
    }
}


/*
 * Compute the side data of primitive variables from characteristic variables.
 */
void
FlowModelBasicUtilitiesSingleSpecies::computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& characteristic_variables,
    const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables)
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
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The number of characteristic variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(characteristic_variables.size()) != d_num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    if (static_cast<int>(projection_variables.size()) != 2)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "There should be two projection variables."
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The interior dimension of the characteristic variables does not match that of patch."
                << std::endl);
        }
    }
    for (int vi = 0; vi < 2; vi++)
    {
        const hier::IntVector interior_dims_projection_var = projection_variables[vi]->getBox().numberCells();
        if (interior_dims_projection_var != interior_dims)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
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
                << ": FlowModelBasicUtilitiesSingleSpecies::"
                << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
                << "The characteristic variables don't have same ghost cell width."
                << std::endl);
        }
    }
    if (num_ghosts_projection_var != projection_variables[1]->getGhostCellWidth())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The projection variables don't have same ghost cell width."
            << std::endl);
    }
    
    if (num_ghosts_projection_var != num_ghosts_primitive_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of primitive"
            << " variables."
            << std::endl);
    }
    if (num_ghosts_projection_var != num_ghosts_characteristic_var)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBasicUtilitiesSingleSpecies::"
            << "computeSideDataOfPrimitiveVariablesFromCharacteristicVariables()\n"
            << "The ghost cell width of the projection variables does not match that of characteristic"
            << " variables."
            << std::endl);
    }
    
    /*
     * Declare containers to store pointers to different data.
     */
    
    std::vector<Real*> V;
    std::vector<Real*> W;
    V.resize(d_num_eqn);
    W.resize(d_num_eqn);
    
    Real* rho_average = nullptr;
    Real* c_average = nullptr;
    
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
        
        rho_average = projection_variables[0]->getPointer(0);
        c_average = projection_variables[1]->getPointer(0);
        
        HAMERS_PRAGMA_SIMD
        for (int i = -num_ghosts_0_characteristic_var;
             i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
             i++)
        {
            // Compute the linear index.
            const int idx_face = i + num_ghosts_0_characteristic_var;
            
            V[0][idx_face] = Real(1)/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                W[1][idx_face] + Real(1)/(c_average[idx_face]*c_average[idx_face])*W[2][idx_face];
            V[1][idx_face] = -Real(1)/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                Real(1)/(rho_average[idx_face]*c_average[idx_face])*W[2][idx_face];
            V[2][idx_face] = W[0][idx_face] + W[2][idx_face];
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
        
        rho_average = projection_variables[0]->getPointer(0);
        c_average = projection_variables[1]->getPointer(0);
        
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -num_ghosts_0_characteristic_var;
                 i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                 i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1);
                
                V[0][idx_face] = Real(1)/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                    W[1][idx_face] + Real(1)/(c_average[idx_face]*c_average[idx_face])*W[3][idx_face];
                V[1][idx_face] = -Real(1)/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                    Real(1)/(rho_average[idx_face]*c_average[idx_face])*W[3][idx_face];
                V[2][idx_face] = W[2][idx_face];
                V[3][idx_face] = W[0][idx_face] + W[3][idx_face];
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
        
        rho_average = projection_variables[0]->getPointer(1);
        c_average = projection_variables[1]->getPointer(1);
        
        for (int j = -num_ghosts_1_characteristic_var;
             j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
             j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index.
                const int idx_face = (i + num_ghosts_0_characteristic_var) +
                    (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var;
                
                V[0][idx_face] = Real(1)/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] + W[1][idx_face] +
                    Real(1)/(c_average[idx_face]*c_average[idx_face])*W[3][idx_face];
                V[1][idx_face] = W[2][idx_face];
                V[2][idx_face] = -Real(1)/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                    Real(1)/(rho_average[idx_face]*c_average[idx_face])*W[3][idx_face];
                V[3][idx_face] = W[0][idx_face] + W[3][idx_face];
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
        
        rho_average = projection_variables[0]->getPointer(0);
        c_average = projection_variables[1]->getPointer(0);
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -num_ghosts_0_characteristic_var;
                     i < interior_dim_0 + 1 + num_ghosts_0_characteristic_var;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1) +
                        (k + num_ghosts_2_characteristic_var)*(ghostcell_dim_0_characteristic_var + 1)*
                            ghostcell_dim_1_characteristic_var;
                    
                    V[0][idx_face] = Real(1)/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        W[1][idx_face] + Real(1)/(c_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[1][idx_face] = -Real(1)/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        Real(1)/(rho_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[2][idx_face] = W[2][idx_face];
                    V[3][idx_face] = W[3][idx_face];
                    V[4][idx_face] = W[0][idx_face] + W[4][idx_face];
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
        
        rho_average = projection_variables[0]->getPointer(1);
        c_average = projection_variables[1]->getPointer(1);
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -num_ghosts_1_characteristic_var;
                 j < interior_dim_1 + 1 + num_ghosts_1_characteristic_var;
                 j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            (ghostcell_dim_1_characteristic_var + 1);
                    
                    V[0][idx_face] = Real(1)/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        W[1][idx_face] + Real(1)/(c_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[1][idx_face] = W[2][idx_face];
                    V[2][idx_face] = -Real(1)/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        Real(1)/(rho_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[3][idx_face] = W[3][idx_face];
                    V[4][idx_face] = W[0][idx_face] + W[4][idx_face];
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
        
        rho_average = projection_variables[0]->getPointer(2);
        c_average = projection_variables[1]->getPointer(2);
        
        for (int k = -num_ghosts_2_characteristic_var;
             k < interior_dim_2 + 1 + num_ghosts_2_characteristic_var;
             k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_face = (i + num_ghosts_0_characteristic_var) +
                        (j + num_ghosts_1_characteristic_var)*ghostcell_dim_0_characteristic_var +
                        (k + num_ghosts_2_characteristic_var)*ghostcell_dim_0_characteristic_var*
                            ghostcell_dim_1_characteristic_var;
                    
                    V[0][idx_face] = Real(1)/(c_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        W[1][idx_face] + Real(1)/(c_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[1][idx_face] = W[2][idx_face];
                    V[2][idx_face] = W[3][idx_face];
                    V[3][idx_face] = -Real(1)/(rho_average[idx_face]*c_average[idx_face])*W[0][idx_face] +
                        Real(1)/(rho_average[idx_face]*c_average[idx_face])*W[4][idx_face];
                    V[4][idx_face] = W[0][idx_face] + W[4][idx_face];
                }
            }
        }
    }
}
