#include "flow/flow_models/single-species/FlowModelRiemannSolverSingleSpecies.hpp"

/*
 * Compute the convective flux from conservative variables.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxFromConservativeVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_minus,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_plus,
    const DIRECTION::TYPE& direction,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type,
    const hier::Box& domain) const
{
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    const int num_eqn = d_flow_model_tmp->getNumberOfEquations();
    
    if (static_cast<int>(conservative_variables_minus.size()) != num_eqn ||
        static_cast<int>(conservative_variables_plus.size()) != num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverSingleSpecies::"
            << "computeConvectiveFluxFromConservativeVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(conservative_variables_minus[ei]);
        TBOX_ASSERT(conservative_variables_plus[ei]);
    }
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = convective_flux->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux->getBox().numberCells() == interior_dims);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(conservative_variables_minus[ei]->getBox().numberCells() == interior_dims);
        TBOX_ASSERT(conservative_variables_plus[ei]->getBox().numberCells() == interior_dims);
    }
    
    TBOX_ASSERT(conservative_variables_minus[0] == conservative_variables_plus[0]);
    for (int ei = 1; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(conservative_variables_minus[ei]->getGhostCellWidth() ==
            conservative_variables_minus[0]->getGhostCellWidth());
        TBOX_ASSERT(conservative_variables_plus[ei]->getGhostCellWidth() ==
            conservative_variables_plus[0]->getGhostCellWidth());
    }
#endif
    
    boost::shared_ptr<pdat::SideData<double> > velocity;
    
    switch (riemann_solver_type)
    {
        case RIEMANN_SOLVER::HLLC:
        {
            switch (direction)
            {
                case DIRECTION::X_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInXDirectionFromPrimitiveVariablesHLLC(
                        convective_flux,
                        velocity,
                        conservative_variables_minus,
                        conservative_variables_plus,
                        domain,
                        false);
                    
                    break;
                }
                case DIRECTION::Y_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC(
                        convective_flux,
                        velocity,
                        conservative_variables_minus,
                        conservative_variables_plus,
                        domain,
                        false);
                    
                    break;
                }
                case DIRECTION::Z_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC(
                        convective_flux,
                        velocity,
                        conservative_variables_minus,
                        conservative_variables_plus,
                        domain,
                        false);
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": Unknown direction requested."
                        << std::endl);
                }
            }
            
            break;
        }
        case RIEMANN_SOLVER::HLLC_HLL:
        {
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
            << ": Unknown Riemann solver requested."
            << std::endl);
        }
    }
}


/*
 * Compute the convective flux from primitive variables.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxFromPrimitiveVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_minus,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_plus,
    const DIRECTION::TYPE& direction,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type,
    const hier::Box& domain) const
{
}


/*
 * Compute the convective flux and velocity from conservative variables.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityFromConservativeVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_minus,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_plus,
    const DIRECTION::TYPE& direction,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type,
    const hier::Box& domain) const
{
    
}


/*
 * Compute the convective flux and velocity from primitive variables.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityFromPrimitiveVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_minus,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_plus,
    const DIRECTION::TYPE& direction,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type,
    const hier::Box& domain) const
{
    
}


/*
 * Compute the convective flux and velocity in the x-direction from conservative variables with
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC_old(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_L,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_R,
    const hier::Box& domain,
    bool compute_velocity) const
{
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    const int num_eqn = d_flow_model_tmp->getNumberOfEquations();
    
    // Get the box that covers the interior of patch.
    const hier::Box interior_box = convective_flux->getBox();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_convective_flux = convective_flux->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_convective_flux =
        convective_flux->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_conservative_variables =
        conservative_variables_L[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_conservative_variables =
        conservative_variables_L[0]->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        const hier::IntVector num_ghosts_min =
            hier::IntVector::min(num_ghosts_convective_flux, num_ghosts_conservative_variables);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(convective_flux->getGhostBox().contains(domain));
        TBOX_ASSERT(conservative_variables_L[0]->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the conservative variables.
     */
    
    std::vector<double*> Q_L;
    std::vector<double*> Q_R;
    Q_L.reserve(num_eqn);
    Q_R.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        Q_L.push_back(conservative_variables_L[ei]->getPointer(0));
        Q_R.push_back(conservative_variables_R[ei]->getPointer(0));
    }
    
    /*
     * Get the equation of state mixing rules and the thermodynamic properties of the species.
     */
    
    const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        d_flow_model_tmp->getEquationOfStateMixingRules();
    
    const int num_thermo_properties = equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<double> thermo_properties;
    std::vector<double*> thermo_properties_ptr;
    std::vector<const double*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_0_conservative_variables = num_ghosts_conservative_variables[0];
        
        /*
         * Get the pointers to the side data of convective flux and conservative variables.
         */
        
        std::vector<double*> F_x;
        F_x.reserve(num_eqn);
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x.push_back(convective_flux->getPointer(0, ei));
        }
        
        std::vector<double*> Q_x_L;
        std::vector<double*> Q_x_R;
        Q_x_L.reserve(num_eqn);
        Q_x_R.reserve(num_eqn);
        for (int ei = 0; ei < num_eqn; ei++)
        {
            Q_x_L.push_back(conservative_variables_L[ei]->getPointer(0, 0));
            Q_x_R.push_back(conservative_variables_R[ei]->getPointer(0, 0));
        }
        
        /*
         * Allocate temporary data.
         */
        
        const int direction_x[1] = {1};
        
        boost::shared_ptr<pdat::SideData<double> > velocity_x_L(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > velocity_x_R(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > internal_energy_x_L(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));

        boost::shared_ptr<pdat::SideData<double> > internal_energy_x_R(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > pressure_x_L(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > pressure_x_R(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > sound_speed_x_L(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > sound_speed_x_R(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > velocity_x_average(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > sound_speed_x_average(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > wave_speed_x_L(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > wave_speed_x_R(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > wave_speed_x_minus(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > wave_speed_x_plus(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > wave_speed_x_star(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > relative_speed_x_star_L(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > relative_speed_x_star_R(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > conservative_variables_x_star_L(
            new pdat::SideData<double>(interior_box, num_eqn, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > conservative_variables_x_star_R(
            new pdat::SideData<double>(interior_box, num_eqn, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > convective_flux_x_L(
            new pdat::SideData<double>(interior_box, num_eqn, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > convective_flux_x_R(
            new pdat::SideData<double>(interior_box, num_eqn, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        /*
         * Get the pointers to the temporary data.
         */
        
        double* u_x_L = velocity_x_L->getPointer(0, 0);
        double* u_x_R = velocity_x_R->getPointer(0, 0);
        
        double* epsilon_x_L = internal_energy_x_L->getPointer(0, 0);
        double* epsilon_x_R = internal_energy_x_R->getPointer(0, 0);
        
        double* p_x_L = pressure_x_L->getPointer(0, 0);
        double* p_x_R = pressure_x_R->getPointer(0, 0);
        
        double* c_x_L = sound_speed_x_L->getPointer(0, 0);
        double* c_x_R = sound_speed_x_R->getPointer(0, 0);
        
        double* u_x_average = velocity_x_average->getPointer(0, 0);
        
        double* c_x_average = sound_speed_x_average->getPointer(0, 0);
        
        double* s_x_L = wave_speed_x_L->getPointer(0, 0);
        double* s_x_R = wave_speed_x_R->getPointer(0, 0);
        
        double* s_x_minus = wave_speed_x_minus->getPointer(0, 0);
        double* s_x_plus = wave_speed_x_plus->getPointer(0, 0);
        
        double* s_x_star = wave_speed_x_star->getPointer(0, 0);
        
        double* Chi_x_star_L = relative_speed_x_star_L->getPointer(0, 0);
        double* Chi_x_star_R = relative_speed_x_star_R->getPointer(0, 0);
        
        std::vector<double*> Q_x_star_L;
        std::vector<double*> Q_x_star_R;
        Q_x_star_L.reserve(num_eqn);
        Q_x_star_R.reserve(num_eqn);
        for (int ei = 0; ei < num_eqn; ei++)
        {
            Q_x_star_L.push_back(conservative_variables_x_star_L->getPointer(0, ei));
            Q_x_star_R.push_back(conservative_variables_x_star_R->getPointer(0, ei));
        }
        
        std::vector<double*> F_x_L;
        std::vector<double*> F_x_R;
        F_x_L.reserve(num_eqn);
        F_x_R.reserve(num_eqn);
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x_L.push_back(convective_flux_x_L->getPointer(0, ei));
            F_x_R.push_back(convective_flux_x_R->getPointer(0, ei));
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            u_x_L[idx] = Q_x_L[1][idx]/Q_x_L[0][idx];
            u_x_R[idx] = Q_x_R[1][idx]/Q_x_R[0][idx];
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            epsilon_x_L[idx] = Q_x_L[2][idx]/Q_x_L[0][idx] -
                double(1)/double(2)*u_x_L[idx]*u_x_L[idx];
            
            epsilon_x_R[idx] = Q_x_R[2][idx]/Q_x_R[0][idx] -
                double(1)/double(2)*u_x_R[idx]*u_x_R[idx];
        }
        
        d_flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_x_L,
                conservative_variables_L[0],
                internal_energy_x_L,
                thermo_properties_const_ptr,
                0,
                domain);
        
        d_flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_x_R,
                conservative_variables_R[0],
                internal_energy_x_R,
                thermo_properties_const_ptr,
                0,
                domain);
        
        d_flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_L,
                conservative_variables_L[0],
                pressure_x_L,
                thermo_properties_const_ptr,
                0,
                domain);
        
        d_flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_R,
                conservative_variables_R[0],
                pressure_x_R,
                thermo_properties_const_ptr,
                0,
                domain);
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            u_x_average[idx] = double(1)/double(2)*(u_x_L[idx] + u_x_R[idx]);
            c_x_average[idx] = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            s_x_L[idx] = fmin(u_x_average[idx] - c_x_average[idx],
                u_x_L[idx] - c_x_L[idx]);
            
            s_x_R[idx] = fmax(u_x_average[idx] + c_x_average[idx],
                u_x_R[idx] + c_x_R[idx]);
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            s_x_minus[idx] = fmin(double(0), s_x_L[idx]);
            s_x_plus[idx] = fmax(double(0), s_x_R[idx]);
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            s_x_star[idx] = (p_x_R[idx] - p_x_L[idx] +
                Q_x_L[1][idx]*(s_x_L[idx] - u_x_L[idx]) -
                Q_x_R[1][idx]*(s_x_R[idx] - u_x_R[idx]))/
                (Q_x_L[0][idx]*(s_x_L[idx] - u_x_L[idx]) -
                 Q_x_R[0][idx]*(s_x_R[idx] - u_x_R[idx]));
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            Chi_x_star_L[idx] = (s_x_L[idx] - u_x_L[idx])/(s_x_L[idx] - s_x_star[idx]);
            Chi_x_star_R[idx] = (s_x_R[idx] - u_x_R[idx])/(s_x_R[idx] - s_x_star[idx]);
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            Q_x_star_L[0][idx] = Chi_x_star_L[idx]*Q_x_L[0][idx];
            Q_x_star_L[1][idx] = Chi_x_star_L[idx]*Q_x_L[0][idx]*s_x_star[idx];
            Q_x_star_L[2][idx] = Chi_x_star_L[idx]*(Q_x_L[2][idx] +
                (s_x_star[idx] - u_x_L[idx])*(Q_x_L[0][idx]*s_x_star[idx] +
                p_x_L[idx]/(s_x_L[idx] - u_x_L[idx])));
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            Q_x_star_R[0][idx] = Chi_x_star_R[idx]*Q_x_R[0][idx];
            Q_x_star_R[1][idx] = Chi_x_star_R[idx]*Q_x_R[0][idx]*s_x_star[idx];
            Q_x_star_R[2][idx] = Chi_x_star_R[idx]*(Q_x_R[2][idx] +
                (s_x_star[idx] - u_x_R[idx])*(Q_x_R[0][idx]*s_x_star[idx] +
                p_x_R[idx]/(s_x_R[idx] - u_x_R[idx])));
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            F_x_L[0][idx] = Q_x_L[1][idx];
            F_x_L[1][idx] = u_x_L[idx]*Q_x_L[1][idx] + p_x_L[idx];
            F_x_L[2][idx] = u_x_L[idx]*(Q_x_L[2][idx] + p_x_L[idx]);
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            F_x_R[0][idx] = Q_x_R[1][idx];
            F_x_R[1][idx] = u_x_R[idx]*Q_x_R[1][idx] + p_x_R[idx];
            F_x_R[2][idx] = u_x_R[idx]*(Q_x_R[2][idx] + p_x_R[idx]);
        }
        
        for (int ei = 0; ei < num_eqn; ei++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
            {
                // Compute the linear index.
                const int idx_convective_flux = i + num_ghosts_0_convective_flux;
                const int idx = i + num_ghosts_0_conservative_variables;
                
                if (s_x_star[idx] > 0)
                {
                    F_x[ei][idx_convective_flux] = F_x_L[ei][idx] +
                        s_x_minus[idx]*(Q_x_star_L[ei][idx] - Q_x_L[ei][idx]);
                }
                else
                {
                    F_x[ei][idx_convective_flux] = F_x_R[ei][idx] +
                        s_x_plus[idx]*(Q_x_star_R[ei][idx] - Q_x_R[ei][idx]);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        
    }
    else if (d_dim == tbox::Dimension(3))
    {
        
    }
}


/*
 * Compute the convective flux and velocity in the x-direction from conservative variables with
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_L,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_R,
    const hier::Box& domain,
    bool compute_velocity) const
{
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    const int num_eqn = d_flow_model_tmp->getNumberOfEquations();
    
    // Get the box that covers the interior of patch.
    const hier::Box interior_box = convective_flux->getBox();
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_convective_flux = convective_flux->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_convective_flux =
        convective_flux->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_conservative_variables =
        conservative_variables_L[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_conservative_variables =
        conservative_variables_L[0]->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        const hier::IntVector num_ghosts_min =
            hier::IntVector::min(num_ghosts_convective_flux, num_ghosts_conservative_variables);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(convective_flux->getGhostBox().contains(domain));
        TBOX_ASSERT(conservative_variables_L[0]->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the conservative variables.
     */
    
    std::vector<double*> Q_L;
    std::vector<double*> Q_R;
    Q_L.reserve(num_eqn);
    Q_R.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        Q_L.push_back(conservative_variables_L[ei]->getPointer(0));
        Q_R.push_back(conservative_variables_R[ei]->getPointer(0));
    }
    
    /*
     * Get the equation of state mixing rules and the thermodynamic properties of the species.
     */
    
    const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules =
        d_flow_model_tmp->getEquationOfStateMixingRules();
    
    const int num_thermo_properties = equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    std::vector<double> thermo_properties;
    std::vector<double*> thermo_properties_ptr;
    std::vector<const double*> thermo_properties_const_ptr;
    
    thermo_properties.resize(num_thermo_properties);
    thermo_properties_ptr.reserve(num_thermo_properties);
    thermo_properties_const_ptr.reserve(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&thermo_properties[ti]);
        thermo_properties_const_ptr.push_back(&thermo_properties[ti]);
    }
    
    equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_0_conservative_variables = num_ghosts_conservative_variables[0];
        
        /*
         * Get the pointers to the side data of convective flux and conservative variables.
         */
        
        std::vector<double*> F_x;
        F_x.reserve(num_eqn);
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x.push_back(convective_flux->getPointer(0, ei));
        }
        
        std::vector<double*> Q_x_L;
        std::vector<double*> Q_x_R;
        Q_x_L.reserve(num_eqn);
        Q_x_R.reserve(num_eqn);
        for (int ei = 0; ei < num_eqn; ei++)
        {
            Q_x_L.push_back(conservative_variables_L[ei]->getPointer(0, 0));
            Q_x_R.push_back(conservative_variables_R[ei]->getPointer(0, 0));
        }
        
        /*
         * Allocate temporary data.
         */
        
        const int direction_x[1] = {1};
        
        boost::shared_ptr<pdat::SideData<double> > internal_energy_x_L(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));

        boost::shared_ptr<pdat::SideData<double> > internal_energy_x_R(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > pressure_x_L(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > pressure_x_R(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > sound_speed_x_L(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        boost::shared_ptr<pdat::SideData<double> > sound_speed_x_R(
            new pdat::SideData<double>(interior_box, 1, num_ghosts_convective_flux,
                hier::IntVector(d_dim, direction_x)));
        
        /*
         * Get the pointers to the temporary data.
         */
        
        double* epsilon_x_L = internal_energy_x_L->getPointer(0, 0);
        double* epsilon_x_R = internal_energy_x_R->getPointer(0, 0);
        
        double* p_x_L = pressure_x_L->getPointer(0, 0);
        double* p_x_R = pressure_x_R->getPointer(0, 0);
        
        double* c_x_L = sound_speed_x_L->getPointer(0, 0);
        double* c_x_R = sound_speed_x_R->getPointer(0, 0);
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            epsilon_x_L[idx] = (Q_x_L[2][idx] -
                double(1)/double(2)*Q_x_L[1][idx]*Q_x_L[1][idx]/Q_x_L[0][idx])/Q_x_L[0][idx];
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear index.
            const int idx = i + num_ghosts_0_conservative_variables;
            
            epsilon_x_R[idx] = (Q_x_R[2][idx] -
                double(1)/double(2)*Q_x_R[1][idx]*Q_x_R[1][idx]/Q_x_R[0][idx])/Q_x_R[0][idx];
        }
        
        d_flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_x_L,
                conservative_variables_L[0],
                internal_energy_x_L,
                thermo_properties_const_ptr,
                0,
                domain);
        
        d_flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computePressure(
                pressure_x_R,
                conservative_variables_R[0],
                internal_energy_x_R,
                thermo_properties_const_ptr,
                0,
                domain);
        
        d_flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_L,
                conservative_variables_L[0],
                pressure_x_L,
                thermo_properties_const_ptr,
                0,
                domain);
        
        d_flow_model_tmp->getEquationOfStateMixingRules()->getEquationOfState()->
            computeSoundSpeed(
                sound_speed_x_R,
                conservative_variables_R[0],
                pressure_x_R,
                thermo_properties_const_ptr,
                0,
                domain);
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
            // Compute the linear indices.
            const int idx_convective_flux = i + num_ghosts_0_convective_flux;
            const int idx = i + num_ghosts_0_conservative_variables;
            
            const double u_x_L = Q_x_L[1][idx]/Q_x_L[0][idx];
            const double u_x_R = Q_x_R[1][idx]/Q_x_R[0][idx];
            
            const double u_x_average = double(1)/double(2)*(u_x_L + u_x_R);
            const double c_x_average = double(1)/double(2)*(c_x_L[idx] + c_x_R[idx]);
            
            const double s_x_L = fmin(u_x_average - c_x_average, u_x_L - c_x_L[idx]);
            const double s_x_R = fmax(u_x_average + c_x_average, u_x_L + c_x_L[idx]);
            
            const double s_x_minus = fmin(double(0), s_x_L);
            const double s_x_plus  = fmax(double(0), s_x_R);
            
            const double s_x_star = (p_x_R[idx] - p_x_L[idx] +
                Q_x_L[1][idx]*(s_x_L - u_x_L) - Q_x_R[1][idx]*(s_x_R - u_x_R))/
                (Q_x_L[0][idx]*(s_x_L - u_x_L) - Q_x_R[0][idx]*(s_x_R - u_x_R));
            
            double Chi_x_star_LR;
            if (s_x_star > 0)
            {
                Chi_x_star_LR = (s_x_L - u_x_L)/(s_x_L - s_x_star);
            }
            else
            {
                Chi_x_star_LR = (s_x_R - u_x_R)/(s_x_R - s_x_star);
            }
            
            double Q_x_star_LR[3];
            if (s_x_star > 0)
            {
                Q_x_star_LR[0] = Chi_x_star_LR*Q_x_L[0][idx];
                Q_x_star_LR[1] = Chi_x_star_LR*Q_x_L[0][idx]*s_x_star;
                Q_x_star_LR[2] = Chi_x_star_LR*(Q_x_L[2][idx] + (s_x_star - u_x_L)*(Q_x_L[0][idx]*s_x_star +
                    p_x_L[idx]/(s_x_L - u_x_L)));
            }
            else
            {
                Q_x_star_LR[0] = Chi_x_star_LR*Q_x_R[0][idx];
                Q_x_star_LR[1] = Chi_x_star_LR*Q_x_R[0][idx]*s_x_star;
                Q_x_star_LR[2] = Chi_x_star_LR*(Q_x_R[2][idx] + (s_x_star - u_x_R)*(Q_x_R[0][idx]*s_x_star +
                    p_x_R[idx]/(s_x_R - u_x_R)));
            }
            
            double F_x_LR[3];
            if (s_x_star > 0)
            {
                F_x_LR[0] = Q_x_L[1][idx];
                F_x_LR[1] = u_x_L*Q_x_L[1][idx] + p_x_L[idx];
                F_x_LR[2] = u_x_L*(Q_x_L[2][idx] + p_x_L[idx]);
            }
            else
            {
                F_x_LR[0] = Q_x_R[1][idx];
                F_x_LR[1] = u_x_R*Q_x_R[1][idx] + p_x_R[idx];
                F_x_LR[2] = u_x_R*(Q_x_R[2][idx] + p_x_R[idx]);
            }
            
            for (int ei = 0; ei < num_eqn; ei++)
            {
                if (s_x_star > 0)
                {
                    F_x[ei][idx_convective_flux] = F_x_LR[ei] + s_x_minus*(Q_x_star_LR[ei] - Q_x_L[ei][idx]);
                }
                else
                {
                    F_x[ei][idx_convective_flux] = F_x_LR[ei] + s_x_plus*(Q_x_star_LR[ei] - Q_x_R[ei][idx]);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        
    }
    else if (d_dim == tbox::Dimension(3))
    {
        
    }
}


/*
 * Compute the convective flux and velocity in the y-direction from conservative variables with
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_B,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_T,
    const hier::Box& domain,
    bool compute_velocity) const
{
    
}


/*
 * Compute the convective flux and velocity in the z-direction from conservative variables with
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_B,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_F,
    const hier::Box& domain,
    bool compute_velocity) const
{
    
}


/*
 * Compute the convective flux and velocity in the x-direction from primitive variables with
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInXDirectionFromPrimitiveVariablesHLLC(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_L,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_R,
    const hier::Box& domain,
    bool compute_velocity) const
{
    
}


/*
 * Compute the convective flux and velocity in the y-direction from primitive variables with
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_B,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_T,
    const hier::Box& domain,
    bool compute_velocity) const
{
    
}


/*
 * Compute the convective flux and velocity in the z-direction from primitive variables with
 * HLLC Riemann solver.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_B,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_F,
    const hier::Box& domain,
    bool compute_velocity) const
{
    
}
