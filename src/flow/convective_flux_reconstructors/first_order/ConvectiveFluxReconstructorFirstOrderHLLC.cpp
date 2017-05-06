#include "flow/convective_flux_reconstructors/first_order/ConvectiveFluxReconstructorFirstOrderHLLC.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

ConvectiveFluxReconstructorFirstOrderHLLC::ConvectiveFluxReconstructorFirstOrderHLLC(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const int& num_species,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db):
        ConvectiveFluxReconstructor(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            num_species,
            flow_model,
            convective_flux_reconstructor_db)
{
    d_num_conv_ghosts = hier::IntVector::getOne(d_dim);
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
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putString("d_shock_capturing_scheme", "FIRST_ORDER_HLLC");
}


/*
 * Compute the convective flux and source due to splitting of convective term on a patch.
 */
void
ConvectiveFluxReconstructorFirstOrderHLLC::computeConvectiveFluxAndSourceOnPatch(
    hier::Patch& patch,
    const boost::shared_ptr<pdat::SideVariable<double> >& variable_convective_flux,
    const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(time);
    NULL_USE(RK_step_number);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // convective ghost cells.
    hier::Box conv_ghost_box = interior_box;
    conv_ghost_box.grow(d_num_conv_ghosts);
    const hier::IntVector conv_ghostcell_dims = conv_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    // Get the side data of convective flux.
    boost::shared_ptr<pdat::SideData<double> > convective_flux(
        BOOST_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(variable_convective_flux, data_context)));
    
    // Get the cell data of source.
    boost::shared_ptr<pdat::CellData<double> > source(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(variable_source, data_context)));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    
    TBOX_ASSERT(source);
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // Allocate temporary patch data.
    boost::shared_ptr<pdat::SideData<double> > velocity_intercell(
        new pdat::SideData<double>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the pointer to the convective flux side data.
         */
        
        std::vector<double*> F_face_x;
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
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables =
            d_flow_model->getGlobalCellDataConservativeVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        num_subghosts_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        
        std::vector<double*> Q;
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
         * Compute the fluxes in the x direction.
         */
        
        // Declare and initialize containers to store the references to the conservative variables.
        std::vector<boost::reference_wrapper<double> > Q_x_L_ref;
        std::vector<boost::reference_wrapper<double> > Q_x_R_ref;
        Q_x_L_ref.reserve(d_num_eqn);
        Q_x_R_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the flux at faces.
        std::vector<boost::reference_wrapper<double> > F_face_x_ref;
        F_face_x_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the velocity at faces.
        std::vector<boost::reference_wrapper<double> > vel_face_x_ref;
        vel_face_x_ref.reserve(d_dim.getValue());
        
        for (int i = 0; i < interior_dims[0] + 1; i++)
        {
            // Compute the linear index.
            const int idx_face_x = i;
            
            // Initialzie container that stores the references to conserevative variables.
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                // Compute the linear indices.
                const int idx_L = i - 1 + num_subghosts_conservative_var[ei][0];
                const int idx_R = i + num_subghosts_conservative_var[ei][0];
                
                Q_x_L_ref.push_back(boost::ref(Q[ei][idx_L]));
                Q_x_R_ref.push_back(boost::ref(Q[ei][idx_R]));
            }
            
            // Initialize container that stores the references to flux at faces.
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_face_x_ref.push_back(boost::ref(F_face_x[ei][idx_face_x]));
            }
            
            // Initialize container that stores the references to the velocity at faces.
            for (int di = 0; di < d_dim.getValue(); di++)
            {
                vel_face_x_ref.push_back(
                    boost::ref(velocity_intercell->getPointer(0, di)[idx_face_x]));
            }
            
            // Apply the Riemann solver.
            d_flow_model->
                computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables(
                    F_face_x_ref,
                    vel_face_x_ref,
                    Q_x_L_ref,
                    Q_x_R_ref,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            
            Q_x_L_ref.clear();
            Q_x_R_ref.clear();
            F_face_x_ref.clear();
            vel_face_x_ref.clear();
            
            // Mulitply fluxes by dt.
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_face_x[ei][idx_face_x] *= dt;
            }
        }
        
        /*
         * Compute the source.
         */
        
        const std::vector<EQN_FORM::TYPE> eqn_form = d_flow_model->getEquationsForm();
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            if (eqn_form[ei] == EQN_FORM::ADVECTIVE)
            {
                double* S = source->getPointer(ei);
                
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute the linear indices. 
                    const int idx_cell_wghost = i + num_subghosts_conservative_var[ei][0];
                    const int idx_cell_nghost = i;
                    const int idx_face_x_L = i;
                    const int idx_face_x_R = i + 1;
                    
                    const double& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                    const double& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                    
                    S[idx_cell_nghost] += dt*Q[ei][idx_cell_wghost]*(u_R - u_L)/dx[0];
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
        
        std::vector<double*> F_face_x;
        std::vector<double*> F_face_y;
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
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables =
            d_flow_model->getGlobalCellDataConservativeVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        num_subghosts_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        
        std::vector<double*> Q;
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
         * Compute the fluxes in the x direction.
         */
        
        // Declare and initialize containers to store the references to the conservative variables.
        std::vector<boost::reference_wrapper<double> > Q_x_L_ref;
        std::vector<boost::reference_wrapper<double> > Q_x_R_ref;
        Q_x_L_ref.reserve(d_num_eqn);
        Q_x_R_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the flux at faces.
        std::vector<boost::reference_wrapper<double> > F_face_x_ref;
        F_face_x_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the velocity at faces.
        std::vector<boost::reference_wrapper<double> > vel_face_x_ref;
        vel_face_x_ref.reserve(d_dim.getValue());
        
        for (int j = 0; j < interior_dims[1]; j++)
        {
            for (int i = 0; i < interior_dims[0] + 1; i++)
            {
                // Compute the linear index.
                const int idx_face_x = i +
                    j*(interior_dims[0] + 1);
                
                // Initialzie container that stores the references to conserevative variables.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    // Compute the linear indices.
                    const int idx_L = (i - 1 + num_subghosts_conservative_var[ei][0]) +
                        (j + num_subghosts_conservative_var[ei][1])*
                            subghostcell_dims_conservative_var[ei][0];
                    
                    const int idx_R = (i + num_subghosts_conservative_var[ei][0]) +
                        (j + num_subghosts_conservative_var[ei][1])*
                            subghostcell_dims_conservative_var[ei][0];
                    
                    Q_x_L_ref.push_back(boost::ref(Q[ei][idx_L]));
                    Q_x_R_ref.push_back(boost::ref(Q[ei][idx_R]));
                }
                
                // Initialize container that stores the references to flux at faces.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_face_x_ref.push_back(boost::ref(F_face_x[ei][idx_face_x]));
                }
                
                // Initialize container that stores the references to the velocity at faces.
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    vel_face_x_ref.push_back(
                        boost::ref(velocity_intercell->getPointer(0, di)[idx_face_x]));
                }
                
                // Apply the Riemann solver.
                d_flow_model->
                    computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables(
                        F_face_x_ref,
                        vel_face_x_ref,
                        Q_x_L_ref,
                        Q_x_R_ref,
                        DIRECTION::X_DIRECTION,
                        RIEMANN_SOLVER::HLLC);
                
                Q_x_L_ref.clear();
                Q_x_R_ref.clear();
                F_face_x_ref.clear();
                vel_face_x_ref.clear();
                
                // Mulitply fluxes by dt.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_face_x[ei][idx_face_x] *= dt;
                }
            }
        }
        
        /*
         * Compute the fluxes in the y direction.
         */
        
        // Declare and initialize containers to store the references to the conservative variables.
        std::vector<boost::reference_wrapper<double> > Q_y_B_ref;
        std::vector<boost::reference_wrapper<double> > Q_y_T_ref;
        Q_y_B_ref.reserve(d_num_eqn);
        Q_y_T_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the flux at faces.
        std::vector<boost::reference_wrapper<double> > F_face_y_ref;
        F_face_y_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the velocity at faces.
        std::vector<boost::reference_wrapper<double> > vel_face_y_ref;
        vel_face_y_ref.reserve(d_dim.getValue());
        
        for (int j = 0; j < interior_dims[1] + 1; j++)
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                // Compute the linear index.
                const int idx_face_y = i +
                    j*interior_dims[0];
                
                // Initialzie container that stores the references to conserevative variables.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    // Compute the linear indices.
                    const int idx_B = (i + num_subghosts_conservative_var[ei][0]) +
                        (j - 1 + num_subghosts_conservative_var[ei][1])*
                            subghostcell_dims_conservative_var[ei][0];
                    
                    const int idx_T = (i + num_subghosts_conservative_var[ei][0]) +
                        (j + num_subghosts_conservative_var[ei][1])*
                            subghostcell_dims_conservative_var[ei][0];
                    
                    Q_y_B_ref.push_back(boost::ref(Q[ei][idx_B]));
                    Q_y_T_ref.push_back(boost::ref(Q[ei][idx_T]));
                }
                
                // Initialize container that stores the references to flux at faces.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_face_y_ref.push_back(boost::ref(F_face_y[ei][idx_face_y]));
                }
                
                // Initialize container that stores the references to the velocity at faces.
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    vel_face_y_ref.push_back(
                        boost::ref(velocity_intercell->getPointer(1, di)[idx_face_y]));
                }
                
                // Apply the Riemann solver.
                d_flow_model->
                    computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables(
                        F_face_y_ref,
                        vel_face_y_ref,
                        Q_y_B_ref,
                        Q_y_T_ref,
                        DIRECTION::Y_DIRECTION,
                        RIEMANN_SOLVER::HLLC);
                
                Q_y_B_ref.clear();
                Q_y_T_ref.clear();
                F_face_y_ref.clear();
                vel_face_y_ref.clear();
                
                // Mulitply fluxes by dt.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_face_y[ei][idx_face_y] *= dt;
                }
            }
        }
        
        /*
         * Compute the source.
         */
        
        const std::vector<EQN_FORM::TYPE> eqn_form = d_flow_model->getEquationsForm();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            if (eqn_form[ei] == EQN_FORM::ADVECTIVE)
            {
                double* S = source->getPointer(ei);
                
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
                        
                        const double& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                        const double& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                        
                        const double& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                        const double& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                        
                        S[idx_cell_nghost] += dt*Q[ei][idx_cell_wghost]*((u_R - u_L)/dx[0] + (v_T - v_B)/dx[1]);
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
        
        std::vector<double*> F_face_x;
        std::vector<double*> F_face_y;
        std::vector<double*> F_face_z;
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
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables =
            d_flow_model->getGlobalCellDataConservativeVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        num_subghosts_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        
        std::vector<double*> Q;
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
         * Compute the fluxes in the x direction.
         */
        
        // Declare and initialize containers to store the references to the conservative variables.
        std::vector<boost::reference_wrapper<double> > Q_x_L_ref;
        std::vector<boost::reference_wrapper<double> > Q_x_R_ref;
        Q_x_L_ref.reserve(d_num_eqn);
        Q_x_R_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the flux at faces.
        std::vector<boost::reference_wrapper<double> > F_face_x_ref;
        F_face_x_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the velocity at faces.
        std::vector<boost::reference_wrapper<double> > vel_face_x_ref;
        vel_face_x_ref.reserve(d_dim.getValue());
        
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
                    
                    // Initialzie container that stores the references to conserevative variables.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        // Compute the linear indices.
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
                        
                        Q_x_L_ref.push_back(boost::ref(Q[ei][idx_L]));
                        Q_x_R_ref.push_back(boost::ref(Q[ei][idx_R]));
                    }
                    
                    // Initialize container that stores the references to flux at faces.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_face_x_ref.push_back(boost::ref(F_face_x[ei][idx_face_x]));
                    }
                    
                    // Initialize container that stores the references to the velocity at faces.
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        vel_face_x_ref.push_back(
                            boost::ref(velocity_intercell->getPointer(0, di)[idx_face_x]));
                    }
                    
                    // Apply the Riemann solver.
                    d_flow_model->
                        computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables(
                            F_face_x_ref,
                            vel_face_x_ref,
                            Q_x_L_ref,
                            Q_x_R_ref,
                            DIRECTION::X_DIRECTION,
                            RIEMANN_SOLVER::HLLC);
                    
                    Q_x_L_ref.clear();
                    Q_x_R_ref.clear();
                    F_face_x_ref.clear();
                    vel_face_x_ref.clear();
                    
                    // Mulitply fluxes by dt.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_face_x[ei][idx_face_x] *= dt;
                    }
                }
            }
        }
        
        /*
         * Compute the fluxes in the y direction.
         */
        
        // Declare and initialize containers to store the references to the conservative variables.
        std::vector<boost::reference_wrapper<double> > Q_y_B_ref;
        std::vector<boost::reference_wrapper<double> > Q_y_T_ref;
        Q_y_B_ref.reserve(d_num_eqn);
        Q_y_T_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the flux at faces.
        std::vector<boost::reference_wrapper<double> > F_face_y_ref;
        F_face_y_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the velocity at faces.
        std::vector<boost::reference_wrapper<double> > vel_face_y_ref;
        vel_face_y_ref.reserve(d_dim.getValue());
        
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
                    
                    // Initialzie container that stores the references to conserevative variables.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
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
                        
                        Q_y_B_ref.push_back(boost::ref(Q[ei][idx_B]));
                        Q_y_T_ref.push_back(boost::ref(Q[ei][idx_T]));
                    }
                    
                    // Initialize container that stores the references to flux at faces.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_face_y_ref.push_back(boost::ref(F_face_y[ei][idx_face_y]));
                    }
                    
                    // Initialize container that stores the references to the velocity at faces.
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        vel_face_y_ref.push_back(
                            boost::ref(velocity_intercell->getPointer(1, di)[idx_face_y]));
                    }
                    
                    // Apply the Riemann solver.
                    d_flow_model->
                        computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables(
                            F_face_y_ref,
                            vel_face_y_ref,
                            Q_y_B_ref,
                            Q_y_T_ref,
                            DIRECTION::Y_DIRECTION,
                            RIEMANN_SOLVER::HLLC);
                    
                    Q_y_B_ref.clear();
                    Q_y_T_ref.clear();
                    F_face_y_ref.clear();
                    vel_face_y_ref.clear();
                    
                    // Mulitply fluxes by dt.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_face_y[ei][idx_face_y] *= dt;
                    }
                }
            }
        }
        
        /*
         * Compute the fluxes in the z direction.
         */
        
        // Declare and initialize containers to store the references to the conservative variables.
        std::vector<boost::reference_wrapper<double> > Q_z_B_ref;
        std::vector<boost::reference_wrapper<double> > Q_z_F_ref;
        Q_z_B_ref.reserve(d_num_eqn);
        Q_z_F_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the flux at faces.
        std::vector<boost::reference_wrapper<double> > F_face_z_ref;
        F_face_z_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the velocity at faces.
        std::vector<boost::reference_wrapper<double> > vel_face_z_ref;
        vel_face_z_ref.reserve(d_dim.getValue());
        
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
                    
                    // Initialzie container that stores the references to conserevative variables.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        // Compute the linear indices.
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
                        
                        Q_z_B_ref.push_back(boost::ref(Q[ei][idx_B]));
                        Q_z_F_ref.push_back(boost::ref(Q[ei][idx_F]));
                    }
                    
                    // Initialize container that stores the references to flux at faces.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_face_z_ref.push_back(boost::ref(F_face_z[ei][idx_face_z]));
                    }
                    
                    // Initialize container that stores the references to the velocity at faces.
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        vel_face_z_ref.push_back(
                            boost::ref(velocity_intercell->getPointer(2, di)[idx_face_z]));
                    }
                    
                    // Apply the Riemann solver.
                    d_flow_model->
                        computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables(
                            F_face_z_ref,
                            vel_face_z_ref,
                            Q_z_B_ref,
                            Q_z_F_ref,
                            DIRECTION::Z_DIRECTION,
                            RIEMANN_SOLVER::HLLC);
                    
                    Q_z_B_ref.clear();
                    Q_z_F_ref.clear();
                    F_face_z_ref.clear();
                    vel_face_z_ref.clear();
                    
                    // Mulitply fluxes by dt.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_face_z[ei][idx_face_z] *= dt;
                    }
                }
            }
        }
        
        /*
         * Compute the source.
         */
        
        const std::vector<EQN_FORM::TYPE> eqn_form = d_flow_model->getEquationsForm();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            if (eqn_form[ei] == EQN_FORM::ADVECTIVE)
            {
                double* S = source->getPointer(ei);
                
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
                            
                            const double& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                            const double& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                            
                            const double& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                            const double& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                            
                            const double& w_B = velocity_intercell->getPointer(2, 2)[idx_face_z_B];
                            const double& w_F = velocity_intercell->getPointer(2, 2)[idx_face_z_F];
                            
                            S[idx_cell_nghost] += dt*Q[ei][idx_cell_wghost]*(
                                (u_R - u_L)/dx[0] + (v_T - v_B)/dx[1] + (w_F - w_B)/dx[2]);
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
