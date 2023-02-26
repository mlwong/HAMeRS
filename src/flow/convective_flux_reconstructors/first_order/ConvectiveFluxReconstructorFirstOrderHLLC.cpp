#include "flow/convective_flux_reconstructors/first_order/ConvectiveFluxReconstructorFirstOrderHLLC.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

ConvectiveFluxReconstructorFirstOrderHLLC::ConvectiveFluxReconstructorFirstOrderHLLC(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const FLOW_MODEL::TYPE& flow_model_type,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db):
        ConvectiveFluxReconstructor(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model_type,
            flow_model,
            convective_flux_reconstructor_db)
{
    d_num_conv_ghosts = hier::IntVector::getOne(d_dim);
    
    d_eqn_form = d_flow_model->getEquationsForm();
    d_has_advective_eqn_form = false;
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
        {
            d_has_advective_eqn_form = true;
        }
    }
}

/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorFirstOrderHLLC::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorFirstOrderHLLC object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorFirstOrderHLLC: this = "
       << (ConvectiveFluxReconstructorFirstOrderHLLC *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorFirstOrderHLLC::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putString("d_shock_capturing_scheme", "FIRST_ORDER_HLLC");
}


/*
 * Compute the convective flux and source due to splitting of convective term on a patch.
 */
void
ConvectiveFluxReconstructorFirstOrderHLLC::computeConvectiveFluxAndSourceOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::SideVariable<Real> >& variable_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_source,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(time);
    NULL_USE(RK_step_number);
    
    d_flow_model->setupRiemannSolver();
    HAMERS_SHARED_PTR<FlowModelRiemannSolver> riemann_solver = d_flow_model->getFlowModelRiemannSolver();
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // convective ghost cells.
    hier::Box conv_ghost_box = interior_box;
    conv_ghost_box.grow(d_num_conv_ghosts);
    const hier::IntVector conv_ghostcell_dims = conv_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    // Get the side data of convective flux.
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux(
        HAMERS_SHARED_PTR_CAST<pdat::SideData<Real>, hier::PatchData>(
            patch.getPatchData(variable_convective_flux, data_context)));
    
    // Get the cell data of source.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > source(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
            patch.getPatchData(variable_source, data_context)));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    
    TBOX_ASSERT(source);
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // Allocate temporary patch data.
    HAMERS_SHARED_PTR<pdat::SideData<Real> > velocity_intercell;
    
    if (d_has_advective_eqn_form)
    {
        velocity_intercell.reset(new pdat::SideData<Real>(
            interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the pointer to the convective flux side data.
         */
        
        std::vector<Real*> F_face_x;
        F_face_x.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_face_x.push_back(convective_flux->getPointer(0, ei));
        }
        
        /*
         * Register the patch and data context.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        /*
         * Get the pointers to the conservative variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_variables =
            d_flow_model->getCellDataOfConservativeVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        num_subghosts_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        
        std::vector<Real*> Q;
        Q.reserve(d_num_eqn);
        
        int count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the conservative variable vector is not in the system of equations, ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                Q.push_back(conservative_variables[vi]->getPointer(di));
                num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
                subghostcell_dims_conservative_var.push_back(conservative_variables[vi]->getGhostBox().numberCells());
                
                count_eqn++;
            }
        }
        
        /*
         * Declare temporary data containers for computing the fluxes at cell edges.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > conservative_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > conservative_variables_plus;
        
        conservative_variables_minus.reserve(d_num_eqn);
        conservative_variables_plus.reserve(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            conservative_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getZero(d_dim)));
            
            conservative_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getZero(d_dim)));
        }
        
        /*
         * Initialize temporary data containers for computing the flux in the x-direction.
         */
        
        std::vector<Real*> Q_minus;
        std::vector<Real*> Q_plus;
        Q_minus.resize(d_num_eqn);
        Q_plus.resize(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q_minus[ei] = conservative_variables_minus[ei]->getPointer(0);
            Q_plus[ei] = conservative_variables_plus[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int i = 0; i < interior_dims[0] + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                const int idx_L = i - 1 + num_subghosts_conservative_var[ei][0];
                const int idx_R = i + num_subghosts_conservative_var[ei][0];
                
                Q_minus[ei][idx_face_x] = Q[ei][idx_L];
                Q_plus[ei][idx_face_x] = Q[ei][idx_R];
            }
        }
        
        /*
         * Compute flux in the x-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromConservativeVariables(
                convective_flux,
                velocity_intercell,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromConservativeVariables(
                convective_flux,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        // Multiply flux by dt.
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int i = 0; i < interior_dims[0] + 1; i++)
            {
                // Compute the linear index.
                const int idx_face_x = i;
                
                F_face_x[ei][idx_face_x] *= Real(dt);
            }
        }
        
        /*
         * Compute the source.
         */
        
        if (d_has_advective_eqn_form)
        {
            for (int ei = 0; ei < d_num_eqn; ei ++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute the linear indices. 
                        const int idx_cell_wghost = i + num_subghosts_conservative_var[ei][0];
                        const int idx_cell_nghost = i;
                        const int idx_face_x_L = i;
                        const int idx_face_x_R = i + 1;
                        
                        const Real& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                        const Real& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                        
                        S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(u_R - u_L)/Real(dx[0]);
                    }
                }
            }
        }
        
        /*
         * Unregister the patch in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the pointers to the convective flux side data.
         */
        
        std::vector<Real*> F_face_x;
        std::vector<Real*> F_face_y;
        F_face_x.reserve(d_num_eqn);
        F_face_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_face_x.push_back(convective_flux->getPointer(0, ei));
            F_face_y.push_back(convective_flux->getPointer(1, ei));
        }
        
        /*
         * Register the patch and data context.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        /*
         * Get the pointers to the conservative variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_variables =
            d_flow_model->getCellDataOfConservativeVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        num_subghosts_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        
        std::vector<Real*> Q;
        Q.reserve(d_num_eqn);
        
        int count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the conservative variable vector is not in the system of equations, ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                Q.push_back(conservative_variables[vi]->getPointer(di));
                num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
                subghostcell_dims_conservative_var.push_back(conservative_variables[vi]->getGhostBox().numberCells());
                
                count_eqn++;
            }
        }
        
        /*
         * Declare temporary data containers for computing the fluxes at cell edges.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > conservative_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > conservative_variables_plus;
        
        conservative_variables_minus.reserve(d_num_eqn);
        conservative_variables_plus.reserve(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            conservative_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getZero(d_dim)));
            
            conservative_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getZero(d_dim)));
        }
        
        std::vector<Real*> Q_minus;
        std::vector<Real*> Q_plus;
        Q_minus.resize(d_num_eqn);
        Q_plus.resize(d_num_eqn);
        
        /*
         * Initialize temporary data containers for computing the flux in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q_minus[ei] = conservative_variables_minus[ei]->getPointer(0);
            Q_plus[ei] = conservative_variables_plus[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0] + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i +
                        j*(interior_dims[0] + 1);
                    
                    const int idx_L = (i - 1 + num_subghosts_conservative_var[ei][0]) +
                        (j + num_subghosts_conservative_var[ei][1])*
                            subghostcell_dims_conservative_var[ei][0];
                    
                    const int idx_R = (i + num_subghosts_conservative_var[ei][0]) +
                        (j + num_subghosts_conservative_var[ei][1])*
                            subghostcell_dims_conservative_var[ei][0];
                    
                    Q_minus[ei][idx_face_x] = Q[ei][idx_L];
                    Q_plus[ei][idx_face_x] = Q[ei][idx_R];
                }
            }
        }
        
        /*
         * Compute flux in the x-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromConservativeVariables(
                convective_flux,
                velocity_intercell,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromConservativeVariables(
                convective_flux,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        // Multiply flux by dt.
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0] + 1; i++)
                {
                    // Compute the linear index.
                    const int idx_face_x = i +
                        j*(interior_dims[0] + 1);
                    
                    F_face_x[ei][idx_face_x] *= Real(dt);
                }
            }
        }
        
        /*
         * Initialize temporary data containers for computing the flux in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q_minus[ei] = conservative_variables_minus[ei]->getPointer(1);
            Q_plus[ei] = conservative_variables_plus[ei]->getPointer(1);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int j = 0; j < interior_dims[1] + 1; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_y = i +
                        j*interior_dims[0];
                    
                    const int idx_B = (i + num_subghosts_conservative_var[ei][0]) +
                        (j - 1 + num_subghosts_conservative_var[ei][1])*
                            subghostcell_dims_conservative_var[ei][0];
                    
                    const int idx_T = (i + num_subghosts_conservative_var[ei][0]) +
                        (j + num_subghosts_conservative_var[ei][1])*
                            subghostcell_dims_conservative_var[ei][0];
                    
                    Q_minus[ei][idx_face_y] = Q[ei][idx_B];
                    Q_plus[ei][idx_face_y] = Q[ei][idx_T];
                }
            }
        }
        
        /*
         * Compute flux in the y-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromConservativeVariables(
                convective_flux,
                velocity_intercell,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromConservativeVariables(
                convective_flux,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        // Multiply flux by dt.
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int j = 0; j < interior_dims[1] + 1; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute the linear index.
                    const int idx_face_y = i +
                        j*interior_dims[0];
                    
                    F_face_y[ei][idx_face_y] *= Real(dt);
                }
            }
        }
        
        /*
         * Compute the source.
         */
        
        if (d_has_advective_eqn_form)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = (i + num_subghosts_conservative_var[ei][0]) +
                                (j + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0];
                            
                            const int idx_cell_nghost = i + j*interior_dims[0];
                            
                            const int idx_face_x_L = i +
                                j*(interior_dims[0] + 1);
                            
                            const int idx_face_x_R = (i + 1) +
                                j*(interior_dims[0] + 1);
                            
                            const int idx_face_y_B = i +
                                j*interior_dims[0];
                            
                            const int idx_face_y_T = i +
                                (j + 1)*interior_dims[0];
                            
                            const Real& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                            const Real& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                            
                            const Real& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                            const Real& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                            
                            S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*((u_R - u_L)/Real(dx[0]) + (v_T - v_B)/Real(dx[1]));
                        }
                    }
                }
            }
        }
        
        /*
         * Unregister the patch in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the pointers to the convective flux side data.
         */
        
        std::vector<Real*> F_face_x;
        std::vector<Real*> F_face_y;
        std::vector<Real*> F_face_z;
        F_face_x.reserve(d_num_eqn);
        F_face_y.reserve(d_num_eqn);
        F_face_z.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_face_x.push_back(convective_flux->getPointer(0, ei));
            F_face_y.push_back(convective_flux->getPointer(1, ei));
            F_face_z.push_back(convective_flux->getPointer(2, ei));
        }
        
        /*
         * Register the patch and data context.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        /*
         * Get the pointers to the conservative variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_variables =
            d_flow_model->getCellDataOfConservativeVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        num_subghosts_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        
        std::vector<Real*> Q;
        Q.reserve(d_num_eqn);
        
        int count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the conservative variable vector is not in the system of equations, ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                Q.push_back(conservative_variables[vi]->getPointer(di));
                num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
                subghostcell_dims_conservative_var.push_back(conservative_variables[vi]->getGhostBox().numberCells());
                
                count_eqn++;
            }
        }
        
        /*
         * Declare temporary data containers for computing the fluxes at cell edges.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > conservative_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > conservative_variables_plus;
        
        conservative_variables_minus.reserve(d_num_eqn);
        conservative_variables_plus.reserve(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            conservative_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getZero(d_dim)));
            
            conservative_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getZero(d_dim)));
        }
        
        std::vector<Real*> Q_minus;
        std::vector<Real*> Q_plus;
        Q_minus.resize(d_num_eqn);
        Q_plus.resize(d_num_eqn);
        
        /*
         * Initialize temporary data containers for computing the flux in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q_minus[ei] = conservative_variables_minus[ei]->getPointer(0);
            Q_plus[ei] = conservative_variables_plus[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int k = 0; k < interior_dims[2]; k++)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0] + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dims[0] + 1) +
                            k*(interior_dims[0] + 1)*interior_dims[1];
                        
                        const int idx_L = (i - 1 + num_subghosts_conservative_var[ei][0]) +
                            (j + num_subghosts_conservative_var[ei][1])*
                                subghostcell_dims_conservative_var[ei][0] +
                            (k + num_subghosts_conservative_var[ei][2])*
                                subghostcell_dims_conservative_var[ei][0]*
                                    subghostcell_dims_conservative_var[ei][1];
                        
                        const int idx_R = (i + num_subghosts_conservative_var[ei][0]) +
                            (j + num_subghosts_conservative_var[ei][1])*
                                subghostcell_dims_conservative_var[ei][0] +
                            (k + num_subghosts_conservative_var[ei][2])*
                                subghostcell_dims_conservative_var[ei][0]*
                                    subghostcell_dims_conservative_var[ei][1];
                        
                        Q_minus[ei][idx_face_x] = Q[ei][idx_L];
                        Q_plus[ei][idx_face_x] = Q[ei][idx_R];
                    }
                }
            }
        }
        
        /*
         * Compute flux in the x-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromConservativeVariables(
                convective_flux,
                velocity_intercell,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromConservativeVariables(
                convective_flux,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        // Multiply flux by dt.
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int k = 0; k < interior_dims[2]; k++)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0] + 1; i++)
                    {
                        // Compute the linear index.
                        const int idx_face_x = i +
                            j*(interior_dims[0] + 1) +
                            k*(interior_dims[0] + 1)*interior_dims[1];
                        
                        F_face_x[ei][idx_face_x] *= Real(dt);
                    }
                }
            }
        }
        
        /*
         * Initialize temporary data containers for computing the flux in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q_minus[ei] = conservative_variables_minus[ei]->getPointer(1);
            Q_plus[ei] = conservative_variables_plus[ei]->getPointer(1);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int k = 0; k < interior_dims[2]; k++)
            {
                for (int j = 0; j < interior_dims[1] + 1; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dims[0] +
                            k*interior_dims[0]*(interior_dims[1] + 1);
                        
                        // Compute the linear indices.
                        const int idx_B = (i + num_subghosts_conservative_var[ei][0]) +
                            (j - 1 + num_subghosts_conservative_var[ei][1])*
                                subghostcell_dims_conservative_var[ei][0] +
                            (k + num_subghosts_conservative_var[ei][2])*
                                subghostcell_dims_conservative_var[ei][0]*
                                    subghostcell_dims_conservative_var[ei][1];
                        
                        const int idx_T = (i + num_subghosts_conservative_var[ei][0]) +
                            (j + num_subghosts_conservative_var[ei][1])*
                                subghostcell_dims_conservative_var[ei][0] +
                            (k + num_subghosts_conservative_var[ei][2])*
                                subghostcell_dims_conservative_var[ei][0]*
                                    subghostcell_dims_conservative_var[ei][1];
                        
                        Q_minus[ei][idx_face_y] = Q[ei][idx_B];
                        Q_plus[ei][idx_face_y] = Q[ei][idx_T];
                    }
                }
            }
        }
        
        /*
         * Compute flux in the y-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromConservativeVariables(
                convective_flux,
                velocity_intercell,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromConservativeVariables(
                convective_flux,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        // Multiply flux by dt.
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int k = 0; k < interior_dims[2]; k++)
            {
                for (int j = 0; j < interior_dims[1] + 1; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute the linear index.
                        const int idx_face_y = i +
                            j*interior_dims[0] +
                            k*interior_dims[0]*(interior_dims[1] + 1);
                        
                        F_face_y[ei][idx_face_y] *= Real(dt);
                    }
                }
            }
        }
        
        /*
         * Initialize temporary data containers for computing the flux in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q_minus[ei] = conservative_variables_minus[ei]->getPointer(2);
            Q_plus[ei] = conservative_variables_plus[ei]->getPointer(2);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int k = 0; k < interior_dims[2] + 1; k++)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_z = i +
                            j*interior_dims[0] +
                            k*interior_dims[0]*interior_dims[1];
                        
                        const int idx_B = (i + num_subghosts_conservative_var[ei][0]) +
                            (j + num_subghosts_conservative_var[ei][1])*
                                subghostcell_dims_conservative_var[ei][0] +
                            (k - 1 + num_subghosts_conservative_var[ei][2])*
                                subghostcell_dims_conservative_var[ei][0]*
                                    subghostcell_dims_conservative_var[ei][1];
                        
                        const int idx_F = (i + num_subghosts_conservative_var[ei][0]) +
                            (j + num_subghosts_conservative_var[ei][1])*
                                subghostcell_dims_conservative_var[ei][0] +
                            (k + num_subghosts_conservative_var[ei][2])*
                                subghostcell_dims_conservative_var[ei][0]*
                                    subghostcell_dims_conservative_var[ei][1];
                        
                        Q_minus[ei][idx_face_z] = Q[ei][idx_B];
                        Q_plus[ei][idx_face_z] = Q[ei][idx_F];
                    }
                }
            }
        }
        
        /*
         * Compute flux in the z-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromConservativeVariables(
                convective_flux,
                velocity_intercell,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::Z_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromConservativeVariables(
                convective_flux,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::Z_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        // Multiply flux by dt.
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int k = 0; k < interior_dims[2] + 1; k++)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute the linear index.
                        const int idx_face_z = i +
                            j*interior_dims[0] +
                            k*interior_dims[0]*interior_dims[1];
                        
                        F_face_z[ei][idx_face_z] *= Real(dt);
                    }
                }
            }
        }
        
        /*
         * Compute the source.
         */
        
        if (d_has_advective_eqn_form)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0]; i++)
                            {
                                // Compute the linear indices. 
                                const int idx_cell_wghost = (i + num_subghosts_conservative_var[ei][0]) +
                                    (j + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0] +
                                    (k + num_subghosts_conservative_var[ei][2])*subghostcell_dims_conservative_var[ei][0]*
                                        subghostcell_dims_conservative_var[ei][1];
                                
                                const int idx_cell_nghost = i +
                                    j*interior_dims[0] +
                                    k*interior_dims[0]*interior_dims[1];
                                
                                const int idx_face_x_L = i +
                                    j*(interior_dims[0] + 1) +
                                    k*(interior_dims[0] + 1)*interior_dims[1];
                                
                                const int idx_face_x_R = (i + 1) +
                                    j*(interior_dims[0] + 1) +
                                    k*(interior_dims[0] + 1)*interior_dims[1];
                                
                                const int idx_face_y_B = i +
                                    j*interior_dims[0] +
                                    k*interior_dims[0]*(interior_dims[1] + 1);
                                
                                const int idx_face_y_T = i +
                                    (j + 1)*interior_dims[0] +
                                    k*interior_dims[0]*(interior_dims[1] + 1);
                                
                                const int idx_face_z_B = i +
                                    j*interior_dims[0] +
                                    k*interior_dims[0]*interior_dims[1];
                                
                                const int idx_face_z_F = i +
                                    j*interior_dims[0] +
                                    (k + 1)*interior_dims[0]*interior_dims[1];
                                
                                const Real& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                                const Real& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                                
                                const Real& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                                const Real& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                                
                                const Real& w_B = velocity_intercell->getPointer(2, 2)[idx_face_z_B];
                                const Real& w_F = velocity_intercell->getPointer(2, 2)[idx_face_z_F];
                                
                                S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                    (u_R - u_L)/Real(dx[0]) + (v_T - v_B)/Real(dx[1]) + (w_F - w_B)/Real(dx[2]));
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Unregister the patch in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(3))
}
