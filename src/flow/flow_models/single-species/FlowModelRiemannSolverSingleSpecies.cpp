#include "flow/flow_models/single-species/FlowModelRiemannSolverSingleSpecies.hpp"

/*
 * Compute the convective flux from conservative variables.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxFromConservativeVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_minus,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_plus,
    const hier::Box& domain,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type)
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
    
    /*
     * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
     */
    
    const hier::IntVector num_ghosts_convective_flux = convective_flux->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_convective_flux =
        convective_flux->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_conservative_variables =
        conservative_variables_minus[0]->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_conservative_variables =
        conservative_variables_minus[0]->getGhostBox().numberCells();
    
    /*
     * Get the local lower indices and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    if (domain.empty())
    {
        hier::IntVector num_ghosts_min(d_dim);
        
        num_ghosts_min = hier::IntVector::min(num_ghosts_convective_flux, num_ghosts_conservative_variables);
        
        hier::Box ghost_box = interior_box;
        ghost_box.grow(num_ghosts_min);
        
        domain_lo = -num_ghosts_min;
        domain_dims = ghost_box.numberCells();
    }
    else
    {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(convective_flux->getGhostBox().contains(domain));
        TBOX_ASSERT(conservative_variables_minus[0]->getGhostBox().contains(domain));
#endif
        
        domain_lo = domain.lower() - interior_box.lower();
        domain_dims = domain.numberCells();
    }
    
    /*
     * Get the pointers to the conservative variables.
     */
    
    std::vector<double*> Q_minus;
    std::vector<double*> Q_plus;
    Q_minus.reserve(num_eqn);
    Q_plus.reserve(num_eqn);
    for (int ei = 0; ei < num_eqn; ei++)
    {
        Q_minus.push_back(conservative_variables_minus[ei]->getPointer(0));
        Q_plus.push_back(conservative_variables_plus[ei]->getPointer(0));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower index, numbers of cells in each dimension and numbers of ghost cells.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        const int num_ghosts_0_convective_flux = num_ghosts_convective_flux[0];
        const int num_ghosts_0_conservative_variables = num_ghosts_conservative_variables[0];
        
        // Get the pointers to the side data of convective flux.
        std::vector<double*> F_x;
        F_x.reserve(num_eqn);
        for (int ei = 0; ei < num_eqn; ei++)
        {
            F_x.push_back(convective_flux->getPointer(0, ei));
        }
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0 + 1; i++)
        {
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
 * Compute the convective flux from primitive variables.
 */
void
FlowModelRiemannSolverSingleSpecies::computeConvectiveFluxFromPrimitiveVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_minus,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_plus,
    const hier::Box& domain,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type)
{
}
