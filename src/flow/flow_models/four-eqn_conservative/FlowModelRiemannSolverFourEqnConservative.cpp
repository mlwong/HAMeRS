#include "flow/flow_models/four-eqn_conservative/FlowModelRiemannSolverFourEqnConservative.hpp"

/*
 * Compute the convective flux from conservative variables.
 */
void
FlowModelRiemannSolverFourEqnConservative::computeConvectiveFluxFromConservativeVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_minus,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_plus,
    const DIRECTION::TYPE& direction,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type,
    const hier::Box& domain) const
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
    
    if (static_cast<int>(conservative_variables_minus.size()) != num_eqn ||
        static_cast<int>(conservative_variables_plus.size()) != num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFourEqnConservative::"
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
    
    TBOX_ASSERT(conservative_variables_minus[0]->getGhostCellWidth() ==
        conservative_variables_plus[0]->getGhostCellWidth());
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
                    computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC(
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
                    computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC(
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
                    computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC(
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
                        << ": FlowModelRiemannSolverFourEqnConservative::"
                        << "computeConvectiveFluxFromConservativeVariables()\n"
                        << "Unknown direction requested."
                        << std::endl);
                }
            }
            
            break;
        }
        case RIEMANN_SOLVER::HLLC_HLL:
        {
            switch (direction)
            {
                case DIRECTION::X_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC_HLL(
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
                    computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC_HLL(
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
                    computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC_HLL(
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
                        << ": FlowModelRiemannSolverFourEqnConservative::"
                        << "computeConvectiveFluxFromConservativeVariables()\n"
                        << "Unknown direction requested."
                        << std::endl);
                }
            }
            
            break;
        } 
        default:
        {
            TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFourEqnConservative::"
            << "computeConvectiveFluxFromConservativeVariables()\n"
            << "Unknown Riemann solver requested."
            << std::endl);
        }
    }
}


/*
 * Compute the convective flux from primitive variables.
 */
void
FlowModelRiemannSolverFourEqnConservative::computeConvectiveFluxFromPrimitiveVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_minus,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_plus,
    const DIRECTION::TYPE& direction,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type,
    const hier::Box& domain) const
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
    
    if (static_cast<int>(primitive_variables_minus.size()) != num_eqn ||
        static_cast<int>(primitive_variables_plus.size()) != num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFourEqnConservative::"
            << "computeConvectiveFluxFromPrimitiveVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(primitive_variables_minus[ei]);
        TBOX_ASSERT(primitive_variables_plus[ei]);
    }
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = convective_flux->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux->getBox().numberCells() == interior_dims);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(primitive_variables_minus[ei]->getBox().numberCells() == interior_dims);
        TBOX_ASSERT(primitive_variables_plus[ei]->getBox().numberCells() == interior_dims);
    }
    
    TBOX_ASSERT(primitive_variables_minus[0]->getGhostCellWidth() ==
        primitive_variables_plus[0]->getGhostCellWidth());
    for (int ei = 1; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(primitive_variables_minus[ei]->getGhostCellWidth() ==
            primitive_variables_minus[0]->getGhostCellWidth());
        TBOX_ASSERT(primitive_variables_plus[ei]->getGhostCellWidth() ==
            primitive_variables_plus[0]->getGhostCellWidth());
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
                        primitive_variables_minus,
                        primitive_variables_plus,
                        domain,
                        false);
                    
                    break;
                }
                case DIRECTION::Y_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC(
                        convective_flux,
                        velocity,
                        primitive_variables_minus,
                        primitive_variables_plus,
                        domain,
                        false);
                    
                    break;
                }
                case DIRECTION::Z_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC(
                        convective_flux,
                        velocity,
                        primitive_variables_minus,
                        primitive_variables_plus,
                        domain,
                        false);
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelRiemannSolverFourEqnConservative::"
                        << "computeConvectiveFluxFromPrimitiveVariables()\n"
                        << "Unknown direction requested."
                        << std::endl);
                }
            }
            
            break;
        }
        case RIEMANN_SOLVER::HLLC_HLL:
        {
            switch (direction)
            {
                case DIRECTION::X_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInXDirectionFromPrimitiveVariablesHLLC_HLL(
                        convective_flux,
                        velocity,
                        primitive_variables_minus,
                        primitive_variables_plus,
                        domain,
                        false);
                    
                    break;
                }
                case DIRECTION::Y_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC_HLL(
                        convective_flux,
                        velocity,
                        primitive_variables_minus,
                        primitive_variables_plus,
                        domain,
                        false);
                    
                    break;
                }
                case DIRECTION::Z_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC_HLL(
                        convective_flux,
                        velocity,
                        primitive_variables_minus,
                        primitive_variables_plus,
                        domain,
                        false);
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelRiemannSolverFourEqnConservative::"
                        << "computeConvectiveFluxFromPrimitiveVariables()\n"
                        << "Unknown direction requested."
                        << std::endl);
                }
            }
            
            break;
        } 
        default:
        {
            TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFourEqnConservative::"
            << "computeConvectiveFluxFromPrimitiveVariables()\n"
            << "Unknown Riemann solver requested."
            << std::endl);
        }
    }
}


/*
 * Compute the convective flux and velocity from conservative variables.
 */
void
FlowModelRiemannSolverFourEqnConservative::computeConvectiveFluxAndVelocityFromConservativeVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_minus,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_plus,
    const DIRECTION::TYPE& direction,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type,
    const hier::Box& domain) const
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
    
    if (static_cast<int>(conservative_variables_minus.size()) != num_eqn ||
        static_cast<int>(conservative_variables_plus.size()) != num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFourEqnConservative::"
            << "computeConvectiveFluxAndVelocityFromConservativeVariables()\n"
            << "The number of conservative variables are incorrect."
            << std::endl);
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    TBOX_ASSERT(velocity);
    
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
    TBOX_ASSERT(velocity->getBox().numberCells() == interior_dims);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(conservative_variables_minus[ei]->getBox().numberCells() == interior_dims);
        TBOX_ASSERT(conservative_variables_plus[ei]->getBox().numberCells() == interior_dims);
    }
    
    TBOX_ASSERT(convective_flux->getGhostCellWidth() ==
        velocity->getGhostCellWidth());
    
    TBOX_ASSERT(conservative_variables_minus[0]->getGhostCellWidth() ==
        conservative_variables_plus[0]->getGhostCellWidth());
    for (int ei = 1; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(conservative_variables_minus[ei]->getGhostCellWidth() ==
            conservative_variables_minus[0]->getGhostCellWidth());
        TBOX_ASSERT(conservative_variables_plus[ei]->getGhostCellWidth() ==
            conservative_variables_plus[0]->getGhostCellWidth());
    }
#endif
    
    switch (riemann_solver_type)
    {
        case RIEMANN_SOLVER::HLLC:
        {
            switch (direction)
            {
                case DIRECTION::X_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC(
                        convective_flux,
                        velocity,
                        conservative_variables_minus,
                        conservative_variables_plus,
                        domain,
                        true);
                    
                    break;
                }
                case DIRECTION::Y_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC(
                        convective_flux,
                        velocity,
                        conservative_variables_minus,
                        conservative_variables_plus,
                        domain,
                        true);
                    
                    break;
                }
                case DIRECTION::Z_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC(
                        convective_flux,
                        velocity,
                        conservative_variables_minus,
                        conservative_variables_plus,
                        domain,
                        true);
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelRiemannSolverFourEqnConservative::"
                        << "computeConvectiveFluxAndVelocityFromConservativeVariables()\n"
                        << "Unknown direction requested."
                        << std::endl);
                }
            }
            
            break;
        }
        case RIEMANN_SOLVER::HLLC_HLL:
        {
            switch (direction)
            {
                case DIRECTION::X_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC_HLL(
                        convective_flux,
                        velocity,
                        conservative_variables_minus,
                        conservative_variables_plus,
                        domain,
                        true);
                    
                    break;
                }
                case DIRECTION::Y_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC_HLL(
                        convective_flux,
                        velocity,
                        conservative_variables_minus,
                        conservative_variables_plus,
                        domain,
                        true);
                    
                    break;
                }
                case DIRECTION::Z_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC_HLL(
                        convective_flux,
                        velocity,
                        conservative_variables_minus,
                        conservative_variables_plus,
                        domain,
                        true);
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelRiemannSolverFourEqnConservative::"
                        << "computeConvectiveFluxAndVelocityFromConservativeVariables()\n"
                        << "Unknown direction requested."
                        << std::endl);
                }
            }
            
            break;
        } 
        default:
        {
            TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFourEqnConservative::"
            << "computeConvectiveFluxAndVelocityFromConservativeVariables()\n"
            << "Unknown Riemann solver requested."
            << std::endl);
        }
    }
}


/*
 * Compute the convective flux and velocity from primitive variables.
 */
void
FlowModelRiemannSolverFourEqnConservative::computeConvectiveFluxAndVelocityFromPrimitiveVariables(
    boost::shared_ptr<pdat::SideData<double> > convective_flux,
    boost::shared_ptr<pdat::SideData<double> > velocity,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_minus,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_plus,
    const DIRECTION::TYPE& direction,
    const RIEMANN_SOLVER::TYPE& riemann_solver_type,
    const hier::Box& domain) const
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
    
    if (static_cast<int>(primitive_variables_minus.size()) != num_eqn ||
        static_cast<int>(primitive_variables_plus.size()) != num_eqn)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFourEqnConservative::"
            << "computeConvectiveFluxAndVelocityFromPrimitiveVariables()\n"
            << "The number of primitive variables are incorrect."
            << std::endl);
    }
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    TBOX_ASSERT(velocity);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(primitive_variables_minus[ei]);
        TBOX_ASSERT(primitive_variables_plus[ei]);
    }
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = convective_flux->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux->getBox().numberCells() == interior_dims);
    TBOX_ASSERT(velocity->getBox().numberCells() == interior_dims);
    
    for (int ei = 0; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(primitive_variables_minus[ei]->getBox().numberCells() == interior_dims);
        TBOX_ASSERT(primitive_variables_plus[ei]->getBox().numberCells() == interior_dims);
    }
    
    TBOX_ASSERT(convective_flux->getGhostCellWidth() ==
        velocity->getGhostCellWidth());
    
    TBOX_ASSERT(primitive_variables_minus[0]->getGhostCellWidth() ==
        primitive_variables_plus[0]->getGhostCellWidth());
    for (int ei = 1; ei < num_eqn; ei++)
    {
        TBOX_ASSERT(primitive_variables_minus[ei]->getGhostCellWidth() ==
            primitive_variables_minus[0]->getGhostCellWidth());
        TBOX_ASSERT(primitive_variables_plus[ei]->getGhostCellWidth() ==
            primitive_variables_plus[0]->getGhostCellWidth());
    }
#endif
    
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
                        primitive_variables_minus,
                        primitive_variables_plus,
                        domain,
                        true);
                    
                    break;
                }
                case DIRECTION::Y_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC(
                        convective_flux,
                        velocity,
                        primitive_variables_minus,
                        primitive_variables_plus,
                        domain,
                        true);
                    
                    break;
                }
                case DIRECTION::Z_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC(
                        convective_flux,
                        velocity,
                        primitive_variables_minus,
                        primitive_variables_plus,
                        domain,
                        true);
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelRiemannSolverFourEqnConservative::"
                        << "computeConvectiveFluxAndVelocityFromPrimitiveVariables()\n"
                        << "Unknown direction requested."
                        << std::endl);
                }
            }
            
            break;
        }
        case RIEMANN_SOLVER::HLLC_HLL:
        {
            switch (direction)
            {
                case DIRECTION::X_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInXDirectionFromPrimitiveVariablesHLLC_HLL(
                        convective_flux,
                        velocity,
                        primitive_variables_minus,
                        primitive_variables_plus,
                        domain,
                        true);
                    
                    break;
                }
                case DIRECTION::Y_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC_HLL(
                        convective_flux,
                        velocity,
                        primitive_variables_minus,
                        primitive_variables_plus,
                        domain,
                        true);
                    
                    break;
                }
                case DIRECTION::Z_DIRECTION:
                {
                    computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC_HLL(
                        convective_flux,
                        velocity,
                        primitive_variables_minus,
                        primitive_variables_plus,
                        domain,
                        true);
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelRiemannSolverFourEqnConservative::"
                        << "computeConvectiveFluxAndVelocityFromPrimitiveVariables()\n"
                        << "Unknown direction requested."
                        << std::endl);
                }
            }
            
            break;
        } 
        default:
        {
            TBOX_ERROR(d_object_name
            << ": FlowModelRiemannSolverFourEqnConservative::"
            << "computeConvectiveFluxAndVelocityFromPrimitiveVariables()\n"
            << "Unknown Riemann solver requested."
            << std::endl);
        }
    }
}
