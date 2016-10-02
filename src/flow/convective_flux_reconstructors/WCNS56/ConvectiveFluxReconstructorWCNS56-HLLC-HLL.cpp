#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS56-HLLC-HLL.hpp"

#define EPSILON 1e-40

ConvectiveFluxReconstructorWCNS56::ConvectiveFluxReconstructorWCNS56(
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
    d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*4;
}


/*
 * Compute the convective fluxes and sources due to hyperbolization
 * of the equations.
 */
void
ConvectiveFluxReconstructorWCNS56::computeConvectiveFluxesAndSources(
    hier::Patch& patch,
    const double time,
    const double dt,
    const int RK_step_number,
    const boost::shared_ptr<pdat::FaceVariable<double> >& variable_convective_flux,
    const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
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
    
    // Get the face data of convective flux.
    boost::shared_ptr<pdat::FaceData<double> > convective_flux(
        BOOST_CAST<pdat::FaceData<double>, hier::PatchData>(
            patch.getPatchData(variable_convective_flux, data_context)));
    
    // Get the cell data of source.
    boost::shared_ptr<pdat::CellData<double> > source(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(variable_source, data_context)));
    
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    
    TBOX_ASSERT(source);
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // Allocate temporary patch data.
    boost::shared_ptr<pdat::FaceData<double> > velocity_intercell(
        new pdat::FaceData<double>(interior_box, d_dim.getValue(), hier::IntVector::getOne(d_dim)));
    
    boost::shared_ptr<pdat::CellData<double> > vorticity_magnitude(
        new pdat::CellData<double>(interior_box, 1, d_num_conv_ghosts));
    
    boost::shared_ptr<pdat::FaceData<double> > convective_flux_midpoint(
        new pdat::FaceData<double>(interior_box, d_num_eqn, hier::IntVector::getOne(d_dim)));

    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRIMITIVE_VARIABLES", d_num_conv_ghosts));
        
        d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
        
        d_flow_model->registerFaceProjectionMatricesOfPrimitiveVariables(
            d_num_conv_ghosts,
            SIMPLE_AVG);
        
        d_flow_model->computeGlobalDerivedCellData();
        
        /*
         * Get the pointers to the velocity and convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        boost::shared_ptr<pdat::CellData<double> > velocity =
            d_flow_model->getGlobalCellData("VELOCITY");
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > convective_flux_node(1);
        convective_flux_node[0] = d_flow_model->getGlobalCellData("CONVECTIVE_FLUX_X");
        
        hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
        hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
        
        double* u = velocity->getPointer(0);
        std::vector<double*> F_x_node;
        F_x_node.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
        }
        std::vector<double*> F_x_midpoint;
        F_x_midpoint.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_x_midpoint.push_back(convective_flux_midpoint->getPointer(0, ei));
        }
        
        /*
         * Get the pointers to the conservative variables and primitive variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables =
            d_flow_model->getGlobalCellDataConservativeVariables();
        std::vector<boost::shared_ptr<pdat::CellData<double> > > primitive_variables =
            d_flow_model->getGlobalCellDataPrimitiveVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        std::vector<hier::IntVector> num_subghosts_primitive_var;
        
        num_subghosts_conservative_var.reserve(d_num_eqn);
        num_subghosts_primitive_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        std::vector<hier::IntVector> subghostcell_dims_primitive_var;
        
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        subghostcell_dims_primitive_var.reserve(d_num_eqn);
        
        std::vector<double*> Q;
        std::vector<double*> V;
        
        Q.reserve(d_num_eqn);
        V.reserve(d_num_eqn);
        
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
        
        count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the primitive variable vector is not in the system of equations, ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                V.push_back(primitive_variables[vi]->getPointer(di));
                num_subghosts_primitive_var.push_back(primitive_variables[vi]->getGhostCellWidth());
                subghostcell_dims_primitive_var.push_back(primitive_variables[vi]->getGhostBox().numberCells());
                
                count_eqn++;
            }
        }
        
        // Declare container to store primitive variables used in WENO interpolation.
        boost::multi_array<const double*, 2> V_array(
            boost::extents[6][d_num_eqn],
            boost::fortran_storage_order());
        
        /*
         * Compute the mid-point fluxes in the x direction.
         */
        
        // Declare indices used in WENO interpolation.
        hier::Index idx_x_L(d_dim);
        hier::Index idx_x_R(d_dim);
        
        // Declare containers to store the WENO interpolated values.
        std::vector<double> V_x_L(d_num_eqn);
        std::vector<double> V_x_R(d_num_eqn);
        
        // Declare and initialize containers to store the references to the
        // WENO interpolate values.
        std::vector<boost::reference_wrapper<double> > V_x_L_ref;
        std::vector<boost::reference_wrapper<double> > V_x_R_ref;
        V_x_L_ref.reserve(d_num_eqn);
        V_x_R_ref.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_x_L_ref.push_back(boost::ref(V_x_L[ei]));
            V_x_R_ref.push_back(boost::ref(V_x_R[ei]));
        }
        
        // Declare container to store the references to the mid-point flux.
        std::vector<boost::reference_wrapper<double> > F_x_midpoint_ref;
        F_x_midpoint_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the mid-point velocity.
        std::vector<boost::reference_wrapper<double> > vel_x_midpoint_ref;
        vel_x_midpoint_ref.reserve(d_dim.getValue());
        
        for (int i = -1; i < interior_dims[0] + 2; i++)
        {
            // Compute the linear index of the face.
            const int idx_face_x = i + 1;
            
            // Compute the indices of left and right cells.
            idx_x_L[0] = i - 1;
            idx_x_R[0] = i;
            
            for (int m = 0; m < 6; m++)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    const int idx_cell = i - 3 + m + num_subghosts_primitive_var[ei][0];
                    V_array[m][ei] = &V[ei][idx_cell];
                }
            }
            
            performWENOInterpolation(
                V_x_L,
                V_x_R,
                V_array,
                idx_x_L,
                idx_x_R,
                X_DIRECTION);
            
            bool is_constant_interpolation = false;
            
            // If the WENO interpolated values are out of bound, use constant interpolation.
            if (!d_flow_model->havePrimitiveVariablesBounded(V_x_L))
            {
                is_constant_interpolation = true;
            }
            
            if (!d_flow_model->havePrimitiveVariablesBounded(V_x_R))
            {
                is_constant_interpolation = true;
            }
            
            if (is_constant_interpolation)
            {
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    // Compute the linear indices of left and right cells.
                    const int idx_x_L = i - 1 + num_subghosts_primitive_var[ei][0];
                    const int idx_x_R = i + num_subghosts_primitive_var[ei][0];
                    
                    V_x_L[ei] = V[ei][idx_x_L];
                    V_x_R[ei] = V[ei][idx_x_R];
                }
            }
            
            // Initialize container that stores the references to the mid-point flux.
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_x_midpoint_ref.push_back(boost::ref(F_x_midpoint[ei][idx_face_x]));
            }
            
            // Initialize container that stores the references to the mid-point velocity.
            for (int di = 0; di < d_dim.getValue(); di++)
            {
                vel_x_midpoint_ref.push_back(
                    boost::ref(velocity_intercell->getPointer(0, di)[idx_face_x]));
            }
            
            // Apply the Riemann solver.
            d_flow_model->
                computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                    F_x_midpoint_ref,
                    vel_x_midpoint_ref,
                    V_x_L_ref,
                    V_x_R_ref,
                    X_DIRECTION,
                    HLLC_RIEMANN_SOLVER);
            
            F_x_midpoint_ref.clear();
            vel_x_midpoint_ref.clear();
        }
        
        /*
         * Compute the fluxes in the x direction.
         */
        
        for (int i = 0; i < interior_dims[0] + 1; i++)
        {
            // Compute the linear indices.
            const int idx_face_x = i;
            const int idx_midpoint_x = i + 1;
            const int idx_node_L = i - 1 + num_subghosts_convective_flux_x[0];
            const int idx_node_R = i + num_subghosts_convective_flux_x[0];
            
            // Compute the fluxes.
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                convective_flux->getPointer(0, ei)[idx_face_x] = dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x + 1] +
                    F_x_midpoint[ei][idx_midpoint_x - 1]) -
                    3.0/10*(F_x_node[ei][idx_node_R] +
                    F_x_node[ei][idx_node_L]) +
                    23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
            }
        }
        
        /*
         * Compute the source.
         */
        
        const std::vector<EQUATION_FORM> eqn_form = d_flow_model->getEquationsForm();
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            if (eqn_form[ei] == ADVECTIVE_EQN)
            {
                double* S = source->getPointer(ei);
                
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute the linear indices. 
                    const int idx_cell_wghost = i + num_subghosts_conservative_var[ei][0];
                    
                    const int idx_cell_wghost_x_L = i - 1 + num_subghosts_velocity[0];
                    
                    const int idx_cell_wghost_x_R = i + 1 + num_subghosts_velocity[0];
                    
                    const int idx_cell_nghost = i;
                    
                    const int idx_face_x_LL = i;
                    
                    const int idx_face_x_L = i + 1;
                    
                    const int idx_face_x_R = i + 2;
                    
                    const int idx_face_x_RR = i + 3;
                    
                    const double& u_LL = velocity_intercell->getPointer(0, 0)[idx_face_x_LL];
                    const double& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                    const double& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                    const double& u_RR = velocity_intercell->getPointer(0, 0)[idx_face_x_RR];
                    
                    S[idx_cell_nghost] += dt*Q[ei][idx_cell_wghost]*(
                        (3.0/2*(u_R - u_L) - 3.0/10*(u[idx_cell_wghost_x_R] - u[idx_cell_wghost_x_L]) +
                         1.0/30*(u_RR - u_LL))/dx[0]);
                }
            }
        }
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("DILATATION", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VORTICITY", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRIMITIVE_VARIABLES", d_num_conv_ghosts));
        
        d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
        
        d_flow_model->registerFaceProjectionMatricesOfPrimitiveVariables(
            d_num_conv_ghosts,
            SIMPLE_AVG);
        
        d_flow_model->computeGlobalDerivedCellData();
        
        /*
         * Get the pointers to the velocity and convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        boost::shared_ptr<pdat::CellData<double> > velocity =
            d_flow_model->getGlobalCellData("VELOCITY");
        
        boost::shared_ptr<pdat::CellData<double> > dilatation =
            d_flow_model->getGlobalCellData("DILATATION");
        
        boost::shared_ptr<pdat::CellData<double> > vorticity =
            d_flow_model->getGlobalCellData("VORTICITY");
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > convective_flux_node(2);
        convective_flux_node[0] = d_flow_model->getGlobalCellData("CONVECTIVE_FLUX_X");
        convective_flux_node[1] = d_flow_model->getGlobalCellData("CONVECTIVE_FLUX_Y");
        
        hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
        hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_dilatation = dilatation->getGhostCellWidth();
        hier::IntVector subghostcell_dims_dilatation = dilatation->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_vorticity = vorticity->getGhostCellWidth();
        hier::IntVector subghostcell_dims_vorticity = vorticity->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
        
        double* u     = velocity->getPointer(0);
        double* v     = velocity->getPointer(1);
        double* theta = dilatation->getPointer(0);
        double* omega = vorticity->getPointer(0);
        double* Omega = vorticity_magnitude->getPointer(0);
        std::vector<double*> F_x_node;
        std::vector<double*> F_y_node;
        F_x_node.reserve(d_num_eqn);
        F_y_node.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
            F_y_node.push_back(convective_flux_node[1]->getPointer(ei));
        }
        std::vector<double*> F_x_midpoint;
        std::vector<double*> F_y_midpoint;
        F_x_midpoint.reserve(d_num_eqn);
        F_y_midpoint.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_x_midpoint.push_back(convective_flux_midpoint->getPointer(0, ei));
            F_y_midpoint.push_back(convective_flux_midpoint->getPointer(1, ei));
        }
        
        // Compute the magnitude of vorticity.
        for (int j = -d_num_conv_ghosts[1] + 1; j < interior_dims[1] + d_num_conv_ghosts[1] - 1; j++)
        {
            for (int i = -d_num_conv_ghosts[0] + 1; i < interior_dims[0] + d_num_conv_ghosts[0] - 1; i++)
            {
                // Compute the linear indices.
                const int idx = (i + d_num_conv_ghosts[0]) +
                    (j + d_num_conv_ghosts[1])*conv_ghostcell_dims[0];
                
                const int idx_vorticity = (i + num_subghosts_vorticity[0]) +
                    (j + num_subghosts_vorticity[1])*subghostcell_dims_vorticity[0];
                
                Omega[idx] = fabs(omega[idx_vorticity]);
            }
        }
        
        /*
         * Get the pointers to the conservative variables and primitive variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables =
            d_flow_model->getGlobalCellDataConservativeVariables();
        std::vector<boost::shared_ptr<pdat::CellData<double> > > primitive_variables =
            d_flow_model->getGlobalCellDataPrimitiveVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        std::vector<hier::IntVector> num_subghosts_primitive_var;
        
        num_subghosts_conservative_var.reserve(d_num_eqn);
        num_subghosts_primitive_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        std::vector<hier::IntVector> subghostcell_dims_primitive_var;
        
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        subghostcell_dims_primitive_var.reserve(d_num_eqn);
        
        std::vector<double*> Q;
        std::vector<double*> V;
        
        Q.reserve(d_num_eqn);
        V.reserve(d_num_eqn);
        
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
        
        count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the primitive variable vector is not in the system of equations, ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                V.push_back(primitive_variables[vi]->getPointer(di));
                num_subghosts_primitive_var.push_back(primitive_variables[vi]->getGhostCellWidth());
                subghostcell_dims_primitive_var.push_back(primitive_variables[vi]->getGhostBox().numberCells());
                
                count_eqn++;
            }
        }
        
        // Declare container to store primitive variables used in WENO interpolation.
        boost::multi_array<const double*, 2> V_array(
            boost::extents[6][d_num_eqn],
            boost::fortran_storage_order());
        
        /*
         * Compute the mid-point fluxes in the x direction.
         */
        
        // Declare indices used in WENO interpolation.
        hier::Index idx_x_L(d_dim);
        hier::Index idx_x_R(d_dim);
        
        // Declare containers to store the WENO interpolated values.
        std::vector<double> V_x_L(d_num_eqn);
        std::vector<double> V_x_R(d_num_eqn);
        
        // Declare and initialize containers to store the references to the
        // WENO interpolate values.
        std::vector<boost::reference_wrapper<double> > V_x_L_ref;
        std::vector<boost::reference_wrapper<double> > V_x_R_ref;
        V_x_L_ref.reserve(d_num_eqn);
        V_x_R_ref.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_x_L_ref.push_back(boost::ref(V_x_L[ei]));
            V_x_R_ref.push_back(boost::ref(V_x_R[ei]));
        }
        
        // Declare container to store the references to the mid-point flux.
        std::vector<boost::reference_wrapper<double> > F_x_midpoint_ref;
        F_x_midpoint_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the mid-point velocity.
        std::vector<boost::reference_wrapper<double> > vel_x_midpoint_ref;
        vel_x_midpoint_ref.reserve(d_dim.getValue());
        
        for (int j = 0; j < interior_dims[1]; j++)
        {
            for (int i = -1; i < interior_dims[0] + 2; i++)
            {
                // Compute the linear index of the face.
                const int idx_face_x = (i + 1) +
                    (j + 1)*(interior_dims[0] + 3);
                
                // Compute the indices of left and right cells.
                idx_x_L[0] = i - 1;
                idx_x_L[1] = j;
                
                idx_x_R[0] = i;
                idx_x_R[1] = j;
                
                for (int m = 0; m < 6; m++)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        const int idx_cell = (i - 3 + m + num_subghosts_primitive_var[ei][0]) +
                            (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0];
                        
                        V_array[m][ei] = &V[ei][idx_cell];
                    }
                }
                
                performWENOInterpolation(
                    V_x_L,
                    V_x_R,
                    V_array,
                    idx_x_L,
                    idx_x_R,
                    X_DIRECTION);
                
                bool is_constant_interpolation = false;
                
                // If the WENO interpolated values are out of bound, use constant interpolation.
                if (!d_flow_model->havePrimitiveVariablesBounded(V_x_L))
                {
                    is_constant_interpolation = true;
                }
                
                if (!d_flow_model->havePrimitiveVariablesBounded(V_x_R))
                {
                    is_constant_interpolation = true;
                }
                
                if (is_constant_interpolation)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        // Compute the linear indices of left and right cells.
                        const int idx_x_L = (i - 1 + num_subghosts_primitive_var[ei][0]) +
                            (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0];
                        
                        const int idx_x_R = (i + num_subghosts_primitive_var[ei][0]) +
                            (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0];
                        
                        V_x_L[ei] = V[ei][idx_x_L];
                        V_x_R[ei] = V[ei][idx_x_R];
                    }
                }
                
                // Compute the average dilatation and magnitude of vorticity.
                const int idx_L = (i - 1 + d_num_conv_ghosts[0]) +
                        (j + d_num_conv_ghosts[1])*conv_ghostcell_dims[0];
                
                const int idx_R = (i + d_num_conv_ghosts[0]) +
                        (j + d_num_conv_ghosts[1])*conv_ghostcell_dims[0];
                
                const int idx_L_dilatation = (i - 1 + num_subghosts_dilatation[0]) +
                        (j + num_subghosts_dilatation[1])*subghostcell_dims_dilatation[0];
                
                const int idx_R_dilatation = (i + num_subghosts_dilatation[0]) +
                        (j + num_subghosts_dilatation[1])*subghostcell_dims_dilatation[0];
                
                const double theta_avg = 0.5*(theta[idx_L_dilatation] + theta[idx_R_dilatation]);
                const double Omega_avg = 0.5*(Omega[idx_L] + Omega[idx_R]);
                
                // Compute the Ducros-like shock sensor.
                const double s = -theta_avg/(fabs(theta_avg) + Omega_avg + EPSILON);
                
                // Initialize container that stores the references to the mid-point flux.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_x_midpoint_ref.push_back(boost::ref(F_x_midpoint[ei][idx_face_x]));
                }
                
                // Initialize container that stores the references to the mid-point velocity.
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    vel_x_midpoint_ref.push_back(
                        boost::ref(velocity_intercell->getPointer(0, di)[idx_face_x]));
                }
                
                // Apply the Riemann solver.
                if (s > 0.65)
                {
                    d_flow_model->
                        computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                            F_x_midpoint_ref,
                            vel_x_midpoint_ref,
                            V_x_L_ref,
                            V_x_R_ref,
                            X_DIRECTION,
                            HLLC_HLL_RIEMANN_SOLVER);
                }
                else
                {
                    d_flow_model->
                        computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                            F_x_midpoint_ref,
                            vel_x_midpoint_ref,
                            V_x_L_ref,
                            V_x_R_ref,
                            X_DIRECTION,
                            HLLC_RIEMANN_SOLVER);
                }
                
                F_x_midpoint_ref.clear();
                vel_x_midpoint_ref.clear();
            }
        }
        
        /*
         * Compute the mid-point fluxes in the y direction.
         */
        
        // Declare indices used in WENO interpolation.
        hier::Index idx_y_B(d_dim);
        hier::Index idx_y_T(d_dim);
        
        // Declare containers to store the WENO interpolated values.
        std::vector<double> V_y_B(d_num_eqn);
        std::vector<double> V_y_T(d_num_eqn);
        
        // Declare and initialize containers to store the references to the
        // WENO interpolate values.
        std::vector<boost::reference_wrapper<double> > V_y_B_ref;
        std::vector<boost::reference_wrapper<double> > V_y_T_ref;
        V_y_B_ref.reserve(d_num_eqn);
        V_y_T_ref.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_y_B_ref.push_back(boost::ref(V_y_B[ei]));
            V_y_T_ref.push_back(boost::ref(V_y_T[ei]));
        }
        
        // Declare container to store the references to the mid-point flux.
        std::vector<boost::reference_wrapper<double> > F_y_midpoint_ref;
        F_y_midpoint_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the mid-point velocity.
        std::vector<boost::reference_wrapper<double> > vel_y_midpoint_ref;
        vel_y_midpoint_ref.reserve(d_dim.getValue());
        
        for (int i = 0; i < interior_dims[0]; i++)
        {
            for (int j = -1; j < interior_dims[1] + 2; j++)
            {
                // Compute the linear index of the face.
                const int idx_face_y = (j + 1) +
                    (i + 1)*(interior_dims[1] + 3);
                
                // Compute the indices of bottom and top cells.
                idx_y_B[0] = i;
                idx_y_B[1] = j - 1;
                
                idx_y_T[0] = i;
                idx_y_T[1] = j;
                
                for (int m = 0; m < 6; m++)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        const int idx_cell = (i + num_subghosts_primitive_var[ei][0]) +
                            (j - 3 + m + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0];
                        
                        V_array[m][ei] = &V[ei][idx_cell];
                    }
                }
                
                performWENOInterpolation(
                    V_y_B,
                    V_y_T,
                    V_array,
                    idx_y_B,
                    idx_y_T,
                    Y_DIRECTION);
                
                bool is_constant_interpolation = false;
                
                // If the WENO interpolated values are out of bound, use constant interpolation.
                if (!d_flow_model->havePrimitiveVariablesBounded(V_y_B))
                {
                    is_constant_interpolation = true;
                }
                
                if (!d_flow_model->havePrimitiveVariablesBounded(V_y_T))
                {
                    is_constant_interpolation = true;
                }
                
                if (is_constant_interpolation)
                {
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        // Compute the indices of bottom and top cells.
                        const int idx_y_B = (i + num_subghosts_primitive_var[ei][0]) +
                            (j - 1 + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0];
                        
                        const int idx_y_T = (i + num_subghosts_primitive_var[ei][0]) +
                            (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0];
                        
                        V_y_B[ei] = V[ei][idx_y_B];
                        V_y_T[ei] = V[ei][idx_y_T];
                    }
                }
                
                // Compute the average dilatation and magnitude of vorticity.
                const int idx_B = (i + d_num_conv_ghosts[0]) +
                        (j - 1 + d_num_conv_ghosts[1])*subghostcell_dims_dilatation[0];
                
                const int idx_T = (i + d_num_conv_ghosts[0]) +
                        (j + d_num_conv_ghosts[1])*subghostcell_dims_dilatation[0];
                
                const int idx_B_dilatation = (i + num_subghosts_dilatation[0]) +
                        (j - 1 + num_subghosts_dilatation[1])*subghostcell_dims_dilatation[0];
                
                const int idx_T_dilatation = (i + num_subghosts_dilatation[0]) +
                        (j + num_subghosts_dilatation[1])*subghostcell_dims_dilatation[0];
                
                const double theta_avg = 0.5*(theta[idx_B_dilatation] + theta[idx_T_dilatation]);
                const double Omega_avg = 0.5*(Omega[idx_B] + Omega[idx_T]);
                
                // Compute the Ducros-like shock sensor.
                const double s = -theta_avg/(fabs(theta_avg) + Omega_avg + EPSILON);
                
                // Initialize container that stores the references to the mid-point flux.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_y_midpoint_ref.push_back(boost::ref(F_y_midpoint[ei][idx_face_y]));
                }
                
                // Initialize container that stores the references to the mid-point velocity.
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    vel_y_midpoint_ref.push_back(
                        boost::ref(velocity_intercell->getPointer(1, di)[idx_face_y]));
                }
                
                // Apply the Riemann solver.
                if (s > 0.65)
                {
                    d_flow_model->
                        computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                            F_y_midpoint_ref,
                            vel_y_midpoint_ref,
                            V_y_B_ref,
                            V_y_T_ref,
                            Y_DIRECTION,
                            HLLC_HLL_RIEMANN_SOLVER);
                }
                else
                {
                    d_flow_model->
                        computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                            F_y_midpoint_ref,
                            vel_y_midpoint_ref,
                            V_y_B_ref,
                            V_y_T_ref,
                            Y_DIRECTION,
                            HLLC_RIEMANN_SOLVER);
                }
                
                F_y_midpoint_ref.clear();
                vel_y_midpoint_ref.clear();
            }
        }
        
        /*
         * Compute the fluxes in the x direction.
         */
        
        for (int j = 0; j < interior_dims[1]; j++)
        {
            for (int i = 0; i < interior_dims[0] + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i +
                    j*(interior_dims[0] + 1);
                
                const int idx_midpoint_x = (i + 1) +
                    (j + 1)*(interior_dims[0] + 3);
                
                const int idx_node_L = (i - 1 + num_subghosts_convective_flux_x[0]) +
                    (j + num_subghosts_convective_flux_x[1])*subghostcell_dims_convective_flux_x[0];
                
                const int idx_node_R = (i + num_subghosts_convective_flux_x[0]) +
                    (j + num_subghosts_convective_flux_x[1])*subghostcell_dims_convective_flux_x[0];
                
                // Compute the fluxes.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    convective_flux->getPointer(0, ei)[idx_face_x] =
                        dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x + 1] +
                        F_x_midpoint[ei][idx_midpoint_x - 1]) -
                        3.0/10*(F_x_node[ei][idx_node_R] +
                        F_x_node[ei][idx_node_L]) +
                        23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
                }
            }
        }
        
        /*
         * Compute the fluxes in the y direction.
         */
        
        for (int i = 0; i < interior_dims[0]; i++)
        {
            for (int j = 0; j < interior_dims[1] + 1; j++)
            {
                // Compute the linear indices.
                const int idx_face_y = j +
                    i*(interior_dims[1] + 1);
                
                const int idx_midpoint_y = (j + 1) +
                    (i + 1)*(interior_dims[1] + 3);
                
                const int idx_node_B = (i + num_subghosts_convective_flux_y[0]) +
                    (j - 1 + num_subghosts_convective_flux_y[1])*subghostcell_dims_convective_flux_y[0];
                
                const int idx_node_T = (i + num_subghosts_convective_flux_y[0]) +
                    (j + num_subghosts_convective_flux_y[1])*subghostcell_dims_convective_flux_y[0];
                
                // Compute the fluxes.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    convective_flux->getPointer(1, ei)[idx_face_y] =
                        dt*(1.0/30*(F_y_midpoint[ei][idx_midpoint_y + 1] +
                        F_y_midpoint[ei][idx_midpoint_y - 1]) -
                        3.0/10*(F_y_node[ei][idx_node_T] +
                        F_y_node[ei][idx_node_B]) +
                        23.0/15*F_y_midpoint[ei][idx_midpoint_y]);
                }
            }
        }
        
        /*
         * Compute the source.
         */
        
        const std::vector<EQUATION_FORM> eqn_form = d_flow_model->getEquationsForm();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            if (eqn_form[ei] == ADVECTIVE_EQN)
            {
                double* S = source->getPointer(ei);
                
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute the linear indices.
                        const int idx_cell_wghost = (i + num_subghosts_conservative_var[ei][0]) +
                            (j + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0];
                        
                        const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_velocity[0]) +
                            (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                        
                        const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_velocity[0]) +
                            (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                        
                        const int idx_cell_wghost_y_B = (i + num_subghosts_velocity[0]) +
                            (j - 1 + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                        
                        const int idx_cell_wghost_y_T = (i + num_subghosts_velocity[0]) +
                            (j + 1 + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                        
                        const int idx_cell_nghost = i + j*interior_dims[0];
                        
                        const int idx_face_x_LL = i + 
                            (j + 1)*(interior_dims[0] + 3);
                        
                        const int idx_face_x_L = (i + 1) +
                            (j + 1)*(interior_dims[0] + 3);
                        
                        const int idx_face_x_R = (i + 2) +
                            (j + 1)*(interior_dims[0] + 3);
                        
                        const int idx_face_x_RR = (i + 3) +
                            (j + 1)*(interior_dims[0] + 3);
                        
                        const int idx_face_y_BB = j +
                            (i + 1)*(interior_dims[1] + 3);
                        
                        const int idx_face_y_B = (j + 1) +
                            (i + 1)*(interior_dims[1] + 3);
                        
                        const int idx_face_y_T = (j + 2) +
                            (i + 1)*(interior_dims[1] + 3);
                        
                        const int idx_face_y_TT = (j + 3) +
                            (i + 1)*(interior_dims[1] + 3);
                        
                        const double& u_LL = velocity_intercell->getPointer(0, 0)[idx_face_x_LL];
                        const double& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                        const double& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                        const double& u_RR = velocity_intercell->getPointer(0, 0)[idx_face_x_RR];
                        
                        const double& v_BB = velocity_intercell->getPointer(1, 1)[idx_face_y_BB];
                        const double& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                        const double& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                        const double& v_TT = velocity_intercell->getPointer(1, 1)[idx_face_y_TT];
                        
                        S[idx_cell_nghost] += dt*Q[ei][idx_cell_wghost]*(
                            (3.0/2*(u_R - u_L) - 3.0/10*(u[idx_cell_wghost_x_R] - u[idx_cell_wghost_x_L]) +
                             1.0/30*(u_RR - u_LL))/dx[0] +
                            (3.0/2*(v_T - v_B) - 3.0/10*(v[idx_cell_wghost_y_T] - v[idx_cell_wghost_y_B]) +
                             1.0/30*(v_TT - v_BB))/dx[1]);
                    }
                }
            }
        }
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("DILATATION", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VORTICITY", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Z", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRIMITIVE_VARIABLES", d_num_conv_ghosts));
        
        d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
        
        d_flow_model->registerFaceProjectionMatricesOfPrimitiveVariables(
            d_num_conv_ghosts,
            SIMPLE_AVG);
        
        d_flow_model->computeGlobalDerivedCellData();
        
        /*
         * Get the pointers to the velocity and convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        boost::shared_ptr<pdat::CellData<double> > velocity =
            d_flow_model->getGlobalCellData("VELOCITY");
        
        boost::shared_ptr<pdat::CellData<double> > dilatation =
            d_flow_model->getGlobalCellData("DILATATION");
        
        boost::shared_ptr<pdat::CellData<double> > vorticity =
            d_flow_model->getGlobalCellData("VORTICITY");
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > convective_flux_node(3);
        convective_flux_node[0] = d_flow_model->getGlobalCellData("CONVECTIVE_FLUX_X");
        convective_flux_node[1] = d_flow_model->getGlobalCellData("CONVECTIVE_FLUX_Y");
        convective_flux_node[2] = d_flow_model->getGlobalCellData("CONVECTIVE_FLUX_Z");
        
        hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
        hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_dilatation = dilatation->getGhostCellWidth();
        hier::IntVector subghostcell_dims_dilatation = dilatation->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_vorticity = vorticity->getGhostCellWidth();
        hier::IntVector subghostcell_dims_vorticity = vorticity->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_z = convective_flux_node[2]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_z = convective_flux_node[2]->getGhostBox().numberCells();
        
        double* u     = velocity->getPointer(0);
        double* v     = velocity->getPointer(1);
        double* w     = velocity->getPointer(2);
        double* theta = dilatation->getPointer(0);
        double* omega_x = vorticity->getPointer(0);
        double* omega_y = vorticity->getPointer(1);
        double* omega_z = vorticity->getPointer(2);
        double* Omega = vorticity_magnitude->getPointer(0);
        std::vector<double*> F_x_node;
        std::vector<double*> F_y_node;
        std::vector<double*> F_z_node;
        F_x_node.reserve(d_num_eqn);
        F_y_node.reserve(d_num_eqn);
        F_z_node.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
            F_y_node.push_back(convective_flux_node[1]->getPointer(ei));
            F_z_node.push_back(convective_flux_node[2]->getPointer(ei));
        }
        std::vector<double*> F_x_midpoint;
        std::vector<double*> F_y_midpoint;
        std::vector<double*> F_z_midpoint;
        F_x_midpoint.reserve(d_num_eqn);
        F_y_midpoint.reserve(d_num_eqn);
        F_z_midpoint.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_x_midpoint.push_back(convective_flux_midpoint->getPointer(0, ei));
            F_y_midpoint.push_back(convective_flux_midpoint->getPointer(1, ei));
            F_z_midpoint.push_back(convective_flux_midpoint->getPointer(2, ei));
        }
        
        // Compute the magnitude of vorticity.
        for (int k = -d_num_conv_ghosts[2] + 1; k < interior_dims[2] + d_num_conv_ghosts[2] - 1; k++)
        {
            for (int j = -d_num_conv_ghosts[1] + 1; j < interior_dims[1] + d_num_conv_ghosts[1] - 1; j++)
            {
                for (int i = -d_num_conv_ghosts[0] + 1; i < interior_dims[0] + d_num_conv_ghosts[0] - 1; i++)
                {
                    // Compute the linear indices
                    const int idx = (i + d_num_conv_ghosts[0]) +
                        (j + d_num_conv_ghosts[1])*conv_ghostcell_dims[0] +
                        (k + d_num_conv_ghosts[2])*conv_ghostcell_dims[0]*conv_ghostcell_dims[1];
                    
                    const int idx_vorticity = (i + num_subghosts_vorticity[0]) +
                        (j + num_subghosts_vorticity[1])*subghostcell_dims_vorticity[0] +
                        (k + num_subghosts_vorticity[2])*subghostcell_dims_vorticity[0]*
                            subghostcell_dims_vorticity[1];
                    
                    Omega[idx] = sqrt(omega_x[idx_vorticity]*omega_x[idx_vorticity] +
                        omega_y[idx_vorticity]*omega_y[idx_vorticity] +
                        omega_z[idx_vorticity]*omega_z[idx_vorticity]);
                }
            }
        }
        
        /*
         * Get the pointers to the conservative variables and primitive variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables =
            d_flow_model->getGlobalCellDataConservativeVariables();
        std::vector<boost::shared_ptr<pdat::CellData<double> > > primitive_variables =
            d_flow_model->getGlobalCellDataPrimitiveVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        std::vector<hier::IntVector> num_subghosts_primitive_var;
        
        num_subghosts_conservative_var.reserve(d_num_eqn);
        num_subghosts_primitive_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        std::vector<hier::IntVector> subghostcell_dims_primitive_var;
        
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        subghostcell_dims_primitive_var.reserve(d_num_eqn);
        
        std::vector<double*> Q;
        std::vector<double*> V;
        
        Q.reserve(d_num_eqn);
        V.reserve(d_num_eqn);
        
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
        
        count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the primitive variable vector is not in the system of equations, ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                V.push_back(primitive_variables[vi]->getPointer(di));
                num_subghosts_primitive_var.push_back(primitive_variables[vi]->getGhostCellWidth());
                subghostcell_dims_primitive_var.push_back(primitive_variables[vi]->getGhostBox().numberCells());
                
                count_eqn++;
            }
        }
        
        // Declare container to store primitive variables used in WENO interpolation.
        boost::multi_array<const double*, 2> V_array(
            boost::extents[6][d_num_eqn],
            boost::fortran_storage_order());
        
        /*
         * Compute the mid-point fluxes in the x direction.
         */
        
        // Declare indices used in WENO interpolation.
        hier::Index idx_x_L(d_dim);
        hier::Index idx_x_R(d_dim);
        
        // Declare containers to store the WENO interpolated values.
        std::vector<double> V_x_L(d_num_eqn);
        std::vector<double> V_x_R(d_num_eqn);
        
        // Declare and initialize containers to store the references to the
        // WENO interpolate values.
        std::vector<boost::reference_wrapper<double> > V_x_L_ref;
        std::vector<boost::reference_wrapper<double> > V_x_R_ref;
        V_x_L_ref.reserve(d_num_eqn);
        V_x_R_ref.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_x_L_ref.push_back(boost::ref(V_x_L[ei]));
            V_x_R_ref.push_back(boost::ref(V_x_R[ei]));
        }
        
        // Declare container to store the references to the mid-point flux.
        std::vector<boost::reference_wrapper<double> > F_x_midpoint_ref;
        F_x_midpoint_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the mid-point velocity.
        std::vector<boost::reference_wrapper<double> > vel_x_midpoint_ref;
        vel_x_midpoint_ref.reserve(d_dim.getValue());
        
        for (int k = 0; k < interior_dims[2]; k++)
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = -1; i < interior_dims[0] + 2; i++)
                {
                    // Compute the linear index of the face.
                    const int idx_face_x = (i + 1) +
                        (j + 1)*(interior_dims[0] + 3) +
                        (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                    
                    // Compute the indices of left and right cells.
                    idx_x_L[0] = i - 1;
                    idx_x_L[1] = j;
                    idx_x_L[2] = k;
                    
                    idx_x_R[0] = i;
                    idx_x_R[1] = j;
                    idx_x_R[2] = k;
                    
                    for (int m = 0; m < 6; m++)
                    {
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            const int idx_cell = (i - 3 + m + num_subghosts_primitive_var[ei][0]) +
                                (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0] +
                                (k + num_subghosts_primitive_var[ei][2])*subghostcell_dims_primitive_var[ei][0]*
                                    subghostcell_dims_primitive_var[ei][1];
                            
                            V_array[m][ei] = &V[ei][idx_cell];
                        }
                    }
                    
                    performWENOInterpolation(
                        V_x_L,
                        V_x_R,
                        V_array,
                        idx_x_L,
                        idx_x_R,
                        X_DIRECTION);
                    
                    bool is_constant_interpolation = false;
                    
                    // If the WENO interpolated values are out of bound, use constant interpolation.
                    if (!d_flow_model->havePrimitiveVariablesBounded(V_x_L))
                    {
                        is_constant_interpolation = true;
                    }
                    
                    if (!d_flow_model->havePrimitiveVariablesBounded(V_x_R))
                    {
                        is_constant_interpolation = true;
                    }
                    
                    if (is_constant_interpolation)
                    {
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            // Compute the linear indices of left cell and right cell.
                            const int idx_x_L = (i - 1 + num_subghosts_primitive_var[ei][0]) +
                                (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0] +
                                (k + num_subghosts_primitive_var[ei][2])*subghostcell_dims_primitive_var[ei][0]*
                                    subghostcell_dims_primitive_var[ei][1];
                            
                            const int idx_x_R = (i + num_subghosts_primitive_var[ei][0]) +
                                (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0] +
                                (k + num_subghosts_primitive_var[ei][2])*subghostcell_dims_primitive_var[ei][0]*
                                    subghostcell_dims_primitive_var[ei][1];
                            
                            V_x_L[ei] = V[ei][idx_x_L];
                            V_x_R[ei] = V[ei][idx_x_R];
                        }
                    }
                    
                    // Compute the average dilatation and magnitude of vorticity.
                    const int idx_L = (i - 1 + d_num_conv_ghosts[0]) +
                        (j + d_num_conv_ghosts[1])*conv_ghostcell_dims[0] +
                        (k + d_num_conv_ghosts[2])*conv_ghostcell_dims[0]*conv_ghostcell_dims[1];
                    
                    const int idx_R = (i + d_num_conv_ghosts[0]) +
                        (j + d_num_conv_ghosts[1])*conv_ghostcell_dims[0] +
                        (k + d_num_conv_ghosts[2])*conv_ghostcell_dims[0]*conv_ghostcell_dims[1];
                    
                    const int idx_L_dilatation = (i - 1 + num_subghosts_dilatation[0]) +
                        (j + num_subghosts_dilatation[1])*subghostcell_dims_dilatation[0] +
                        (k + num_subghosts_dilatation[2])*subghostcell_dims_dilatation[0]*
                            subghostcell_dims_dilatation[1];
                    
                    const int idx_R_dilatation = (i + num_subghosts_dilatation[0]) +
                        (j + num_subghosts_dilatation[1])*subghostcell_dims_dilatation[0] +
                        (k + num_subghosts_dilatation[2])*subghostcell_dims_dilatation[0]*
                            subghostcell_dims_dilatation[1];
                    
                    const double theta_avg = 0.5*(theta[idx_L_dilatation] + theta[idx_R_dilatation]);
                    const double Omega_avg = 0.5*(Omega[idx_L] + Omega[idx_R]);
                    
                    // Compute the Ducros-like shock sensor.
                    const double s = -theta_avg/(fabs(theta_avg) + Omega_avg + EPSILON);
                    
                    // Initialize container that stores the references to the mid-point flux.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_midpoint_ref.push_back(boost::ref(F_x_midpoint[ei][idx_face_x]));
                    }
                    
                    // Initialize container that stores the references to the mid-point velocity.
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        vel_x_midpoint_ref.push_back(
                            boost::ref(velocity_intercell->getPointer(0, di)[idx_face_x]));
                    }
                    
                    // Apply the Riemann solver.
                    if (s > 0.65)
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_x_midpoint_ref,
                                vel_x_midpoint_ref,
                                V_x_L_ref,
                                V_x_R_ref,
                                X_DIRECTION,
                                HLLC_HLL_RIEMANN_SOLVER);
                    }
                    else
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_x_midpoint_ref,
                                vel_x_midpoint_ref,
                                V_x_L_ref,
                                V_x_R_ref,
                                X_DIRECTION,
                                HLLC_RIEMANN_SOLVER);
                    }
                    
                    F_x_midpoint_ref.clear();
                    vel_x_midpoint_ref.clear();
                }
            }
        }
        
        /*
         * Compute the mid-point fluxes in the y direction.
         */
        
        // Declare indices used in WENO interpolation.
        hier::Index idx_y_B(d_dim);
        hier::Index idx_y_T(d_dim);
        
        // Declare containers to store the WENO interpolated values.
        std::vector<double> V_y_B(d_num_eqn);
        std::vector<double> V_y_T(d_num_eqn);
        
        // Declare and initialize containers to store the references to the
        // WENO interpolate values.
        std::vector<boost::reference_wrapper<double> > V_y_B_ref;
        std::vector<boost::reference_wrapper<double> > V_y_T_ref;
        V_y_B_ref.reserve(d_num_eqn);
        V_y_T_ref.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_y_B_ref.push_back(boost::ref(V_y_B[ei]));
            V_y_T_ref.push_back(boost::ref(V_y_T[ei]));
        }
        
        // Declare container to store the references to the mid-point flux.
        std::vector<boost::reference_wrapper<double> > F_y_midpoint_ref;
        F_y_midpoint_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the mid-point velocity.
        std::vector<boost::reference_wrapper<double> > vel_y_midpoint_ref;
        vel_y_midpoint_ref.reserve(d_dim.getValue());
        
        for (int i = 0; i < interior_dims[0]; i++)
        {
            for (int k = 0; k < interior_dims[2]; k++)
            {
                for (int j = -1; j < interior_dims[1] + 2; j++)
                {
                    // Compute the linear index of the face.
                    const int idx_face_y = (j + 1) +
                        (k + 1)*(interior_dims[1] + 3) +
                        (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                    
                    // Compute the indices of bottom and top cells.
                    idx_y_B[0] = i;
                    idx_y_B[1] = j - 1;
                    idx_y_B[2] = k;
                    
                    idx_y_T[0] = i;
                    idx_y_T[1] = j;
                    idx_y_T[2] = k;
                    
                    for (int m = 0; m < 6; m++)
                    {
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            const int idx_cell = (i + num_subghosts_primitive_var[ei][0]) +
                                (j - 3 + m + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0] +
                                (k + num_subghosts_primitive_var[ei][2])*subghostcell_dims_primitive_var[ei][0]*
                                    subghostcell_dims_primitive_var[ei][1];
                            
                            V_array[m][ei] = &V[ei][idx_cell];
                        }
                    }
                    
                    performWENOInterpolation(
                        V_y_B,
                        V_y_T,
                        V_array,
                        idx_y_B,
                        idx_y_T,
                        Y_DIRECTION);
                    
                    bool is_constant_interpolation = false;
                    
                    // If the WENO interpolated values are out of bound, use constant interpolation.
                    if (!d_flow_model->havePrimitiveVariablesBounded(V_y_B))
                    {
                        is_constant_interpolation = true;
                    }
                    
                    if (!d_flow_model->havePrimitiveVariablesBounded(V_y_T))
                    {
                        is_constant_interpolation = true;
                    }
                    
                    if (is_constant_interpolation)
                    {
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            // Compute the indices of bottom cell and top cell.
                            const int idx_y_B = (i + num_subghosts_primitive_var[ei][0]) +
                                (j - 1 + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0] +
                                (k + num_subghosts_primitive_var[ei][2])*subghostcell_dims_primitive_var[ei][0]*
                                    subghostcell_dims_primitive_var[ei][1];
                            
                            const int idx_y_T = (i + num_subghosts_primitive_var[ei][0]) +
                                (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0] +
                                (k + num_subghosts_primitive_var[ei][2])*subghostcell_dims_primitive_var[ei][0]*
                                    subghostcell_dims_primitive_var[ei][1];
                            
                            V_y_B[ei] = V[ei][idx_y_B];
                            V_y_T[ei] = V[ei][idx_y_T];
                        }
                    }
                    
                    // Compute the average dilatation and magnitude of vorticity.
                    const int idx_B = (i + d_num_conv_ghosts[0]) +
                        (j - 1 + d_num_conv_ghosts[1])*conv_ghostcell_dims[0] +
                        (k + d_num_conv_ghosts[2])*conv_ghostcell_dims[0]*conv_ghostcell_dims[1];
                    
                    const int idx_T = (i + d_num_conv_ghosts[0]) +
                        (j + d_num_conv_ghosts[1])*conv_ghostcell_dims[0] +
                        (k + d_num_conv_ghosts[2])*conv_ghostcell_dims[0]*conv_ghostcell_dims[1];
                    
                    const int idx_B_dilatation = (i + num_subghosts_dilatation[0]) +
                        (j - 1 + num_subghosts_dilatation[1])*subghostcell_dims_dilatation[0] +
                        (k + num_subghosts_dilatation[2])*subghostcell_dims_dilatation[0]*
                            subghostcell_dims_dilatation[1];
                    
                    const int idx_T_dilatation = (i + num_subghosts_dilatation[0]) +
                        (j + num_subghosts_dilatation[1])*subghostcell_dims_dilatation[0] +
                        (k + num_subghosts_dilatation[2])*subghostcell_dims_dilatation[0]*
                            subghostcell_dims_dilatation[1];
                    
                    const double theta_avg = 0.5*(theta[idx_B_dilatation] + theta[idx_T_dilatation]);
                    const double Omega_avg = 0.5*(Omega[idx_B] + Omega[idx_T]);
                    
                    // Compute the Ducros-like shock sensor.
                    const double s = -theta_avg/(fabs(theta_avg) + Omega_avg + EPSILON);
                    
                    // Initialize container that stores the references to the mid-point flux.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_y_midpoint_ref.push_back(boost::ref(F_y_midpoint[ei][idx_face_y]));
                    }
                    
                    // Initialize container that stores the references to the mid-point velocity.
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        vel_y_midpoint_ref.push_back(
                            boost::ref(velocity_intercell->getPointer(1, di)[idx_face_y]));
                    }
                    
                    // Apply the Riemann solver.
                    if (s > 0.65)
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_y_midpoint_ref,
                                vel_y_midpoint_ref,
                                V_y_B_ref,
                                V_y_T_ref,
                                Y_DIRECTION,
                                HLLC_HLL_RIEMANN_SOLVER);
                    }
                    else
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_y_midpoint_ref,
                                vel_y_midpoint_ref,
                                V_y_B_ref,
                                V_y_T_ref,
                                Y_DIRECTION,
                                HLLC_RIEMANN_SOLVER);
                    }
                    
                    F_y_midpoint_ref.clear();
                    vel_y_midpoint_ref.clear();
                }
            }
        }
        
        /*
         * Compute the mid-point fluxes in the z direction.
         */
        
        // Declare indices used in WENO interpolation.
        hier::Index idx_z_B(d_dim);
        hier::Index idx_z_F(d_dim);
        
        // Declare containers to store the WENO interpolated values.
        std::vector<double> V_z_B(d_num_eqn);
        std::vector<double> V_z_F(d_num_eqn);
        
        // Declare and initialize containers to store the references to the
        // WENO interpolate values.
        std::vector<boost::reference_wrapper<double> > V_z_B_ref;
        std::vector<boost::reference_wrapper<double> > V_z_F_ref;
        V_z_B_ref.reserve(d_num_eqn);
        V_z_F_ref.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_z_B_ref.push_back(boost::ref(V_z_B[ei]));
            V_z_F_ref.push_back(boost::ref(V_z_F[ei]));
        }
        
        // Declare container to store the references to the mid-point flux.
        std::vector<boost::reference_wrapper<double> > F_z_midpoint_ref;
        F_z_midpoint_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the mid-point velocity.
        std::vector<boost::reference_wrapper<double> > vel_z_midpoint_ref;
        vel_z_midpoint_ref.reserve(d_dim.getValue());
        
        for (int j = 0; j < interior_dims[1]; j++)
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                for (int k = -1; k < interior_dims[2] + 2; k++)
                {
                    // Compute the linear index of the face.
                    const int idx_face_z = (k + 1) +
                        (i + 1)*(interior_dims[2] + 3) +
                        (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                    
                    // Compute the indices of back and front cells.
                    idx_z_B[0] = i;
                    idx_z_B[1] = j;
                    idx_z_B[2] = k - 1;
                    
                    idx_z_F[0] = i;
                    idx_z_F[1] = j;
                    idx_z_F[2] = k;
                    
                    for (int m = 0; m < 6; m++)
                    {
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            const int idx_cell = (i + num_subghosts_primitive_var[ei][0]) +
                                (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0] +
                                (k - 3 + m + num_subghosts_primitive_var[ei][2])*subghostcell_dims_primitive_var[ei][0]*
                                    subghostcell_dims_primitive_var[ei][1];
                            
                            V_array[m][ei] = &V[ei][idx_cell];
                        }
                    }
                    
                    performWENOInterpolation(
                        V_z_B,
                        V_z_F,
                        V_array,
                        idx_z_B,
                        idx_z_F,
                        Z_DIRECTION);
                    
                    bool is_constant_interpolation = false;
                    
                    // If the WENO interpolated values are out of bound, use constant interpolation.
                    if (!d_flow_model->havePrimitiveVariablesBounded(V_z_B))
                    {
                        is_constant_interpolation = true;
                    }
                    
                    if (!d_flow_model->havePrimitiveVariablesBounded(V_z_F))
                    {
                        is_constant_interpolation = true;
                    }
                    
                    if (is_constant_interpolation)
                    {
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            // Compute the indices of back and front cells.
                            const int idx_z_B = (i + num_subghosts_primitive_var[ei][0]) +
                                (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0] +
                                (k - 1 + num_subghosts_primitive_var[ei][2])*subghostcell_dims_primitive_var[ei][0]*
                                    subghostcell_dims_primitive_var[ei][1];
                            
                            const int idx_z_F = (i + num_subghosts_primitive_var[ei][0]) +
                                (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0] +
                                (k + num_subghosts_primitive_var[ei][2])*subghostcell_dims_primitive_var[ei][0]*
                                    subghostcell_dims_primitive_var[ei][1];
                            
                            V_z_B[ei] = V[ei][idx_z_B];
                            V_z_F[ei] = V[ei][idx_z_F];
                        }
                    }
                    
                    // Compute the average dilatation and magnitude of vorticity.
                    const int idx_B = (i + d_num_conv_ghosts[0]) +
                        (j + d_num_conv_ghosts[1])*conv_ghostcell_dims[0] +
                        (k - 1 + d_num_conv_ghosts[2])*conv_ghostcell_dims[0]*conv_ghostcell_dims[1];
                    
                    const int idx_F = (i + d_num_conv_ghosts[0]) +
                        (j + d_num_conv_ghosts[1])*conv_ghostcell_dims[0] +
                        (k + d_num_conv_ghosts[2])*conv_ghostcell_dims[0]*conv_ghostcell_dims[1];
                    
                    const int idx_B_dilatation = (i + num_subghosts_dilatation[0]) +
                        (j + num_subghosts_dilatation[1])*subghostcell_dims_dilatation[0] +
                        (k - 1 + num_subghosts_dilatation[2])*subghostcell_dims_dilatation[0]*
                            subghostcell_dims_dilatation[1];
                    
                    const int idx_F_dilatation = (i + num_subghosts_dilatation[0]) +
                        (j + num_subghosts_dilatation[1])*subghostcell_dims_dilatation[0] +
                        (k + num_subghosts_dilatation[2])*subghostcell_dims_dilatation[0]*
                            subghostcell_dims_dilatation[1];
                    
                    const double theta_avg = 0.5*(theta[idx_B_dilatation] + theta[idx_F_dilatation]);
                    const double Omega_avg = 0.5*(Omega[idx_B] + Omega[idx_F]);
                    
                    // Compute the Ducros-like shock sensor.
                    const double s = -theta_avg/(fabs(theta_avg) + Omega_avg + EPSILON);
                    
                    // Initialize container that stores the references to the mid-point flux.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_z_midpoint_ref.push_back(boost::ref(F_z_midpoint[ei][idx_face_z]));
                    }
                    
                    // Initialize container that stores the references to the mid-point velocity.
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        vel_z_midpoint_ref.push_back(
                            boost::ref(velocity_intercell->getPointer(2, di)[idx_face_z]));
                    }
                    
                    // Apply the Riemann solver.
                    if (s > 0.65)
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_z_midpoint_ref,
                                vel_z_midpoint_ref,
                                V_z_B_ref,
                                V_z_F_ref,
                                Z_DIRECTION,
                                HLLC_HLL_RIEMANN_SOLVER);
                    }
                    else
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_z_midpoint_ref,
                                vel_z_midpoint_ref,
                                V_z_B_ref,
                                V_z_F_ref,
                                Z_DIRECTION,
                                HLLC_RIEMANN_SOLVER);
                    }
                    
                    F_z_midpoint_ref.clear();
                    vel_z_midpoint_ref.clear();
                }
            }
        }
        
        /*
         * Compute the fluxes in the x direction.
         */
        for (int k = 0; k < interior_dims[2]; k++)
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0] + 1; i++)
                {
                    // Compute the indices.
                    const int idx_face_x = i +
                        j*(interior_dims[0] + 1) +
                        k*(interior_dims[0] + 1)*interior_dims[1];
                    
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dims[0] + 3) +
                        (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                    
                    const int idx_node_L = (i - 1 + num_subghosts_convective_flux_x[0]) +
                        (j + num_subghosts_convective_flux_x[1])*subghostcell_dims_convective_flux_x[0] +
                        (k + num_subghosts_convective_flux_x[2])*subghostcell_dims_convective_flux_x[0]*
                            subghostcell_dims_convective_flux_x[1];
                    
                    const int idx_node_R = (i + num_subghosts_convective_flux_x[0]) +
                        (j + num_subghosts_convective_flux_x[1])*subghostcell_dims_convective_flux_x[0] +
                        (k + num_subghosts_convective_flux_x[2])*subghostcell_dims_convective_flux_x[0]*
                            subghostcell_dims_convective_flux_x[1];
                    
                    // Compute the fluxes.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        convective_flux->getPointer(0, ei)[idx_face_x] = dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x + 1] +
                            F_x_midpoint[ei][idx_midpoint_x - 1]) -
                            3.0/10*(F_x_node[ei][idx_node_R] +
                            F_x_node[ei][idx_node_L]) +
                            23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
                    }
                }
            }
        }
        
        /*
         * Compute the fluxes in the y direction.
         */
        
        for (int i = 0; i < interior_dims[0]; i++)
        {
            for (int k = 0; k < interior_dims[2]; k++)
            {
                for (int j = 0; j < interior_dims[1] + 1; j++)
                {
                    // Compute the indices.
                    const int idx_face_y = j +
                        k*(interior_dims[1] + 1) +
                        i*(interior_dims[1] + 1)*interior_dims[2];
                    
                    const int idx_midpoint_y = (j + 1) +
                        (k + 1)*(interior_dims[1] + 3) +
                        (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                    
                    const int idx_node_B = (i + num_subghosts_convective_flux_y[0]) +
                        (j - 1 + num_subghosts_convective_flux_y[1])*subghostcell_dims_convective_flux_y[0] +
                        (k + num_subghosts_convective_flux_y[2])*subghostcell_dims_convective_flux_y[0]*
                            subghostcell_dims_convective_flux_y[1];
                    
                    const int idx_node_T = (i + num_subghosts_convective_flux_y[0]) +
                        (j + num_subghosts_convective_flux_y[1])*subghostcell_dims_convective_flux_y[0] +
                        (k + num_subghosts_convective_flux_y[2])*subghostcell_dims_convective_flux_y[0]*
                            subghostcell_dims_convective_flux_y[1];
                    
                    // Compute the fluxes.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        convective_flux->getPointer(1, ei)[idx_face_y] = 
                            dt*(1.0/30*(F_y_midpoint[ei][idx_midpoint_y + 1] +
                            F_y_midpoint[ei][idx_midpoint_y - 1]) -
                            3.0/10*(F_y_node[ei][idx_node_T] +
                            F_y_node[ei][idx_node_B]) +
                            23.0/15*F_y_midpoint[ei][idx_midpoint_y]);
                    }
                }
            }
        }
        
        /*
         * Compute the fluxes in the z direction.
         */
        
        for (int j = 0; j < interior_dims[1]; j++)
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                for (int k = 0; k < interior_dims[2] + 1; k++)
                {
                    // Compute the indices.
                    const int idx_face_z = k +
                        i*(interior_dims[2] + 1) +
                        j*(interior_dims[2] + 1)*interior_dims[0];
                    
                    const int idx_midpoint_z = (k + 1) +
                        (i + 1)*(interior_dims[2] + 3) +
                        (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                    
                    const int idx_node_B = (i + num_subghosts_convective_flux_z[0]) +
                        (j + num_subghosts_convective_flux_z[1])*subghostcell_dims_convective_flux_z[0] +
                        (k - 1 + num_subghosts_convective_flux_z[2])*subghostcell_dims_convective_flux_z[0]*
                            subghostcell_dims_convective_flux_z[1];
                    
                    const int idx_node_F = (i + num_subghosts_convective_flux_z[0]) +
                        (j + num_subghosts_convective_flux_z[1])*subghostcell_dims_convective_flux_z[0] +
                        (k + num_subghosts_convective_flux_z[2])*subghostcell_dims_convective_flux_z[0]*
                            subghostcell_dims_convective_flux_z[1];
                    
                    // Compute the fluxes.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        convective_flux->getPointer(2, ei)[idx_face_z] = 
                            dt*(1.0/30*(F_z_midpoint[ei][idx_midpoint_z + 1] +
                            F_z_midpoint[ei][idx_midpoint_z - 1]) -
                            3.0/10*(F_z_node[ei][idx_node_F] +
                            F_z_node[ei][idx_node_B]) +
                            23.0/15*F_z_midpoint[ei][idx_midpoint_z]);
                    }
                }
            }
        }
        
        /*
         * Compute the source.
         */
        
        const std::vector<EQUATION_FORM> eqn_form = d_flow_model->getEquationsForm();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            if (eqn_form[ei] == ADVECTIVE_EQN)
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
                            
                            const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_velocity[0]) +
                                (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                                (k + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                    subghostcell_dims_velocity[1];
                            
                            const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_velocity[0]) +
                                (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                                (k + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                    subghostcell_dims_velocity[1];
                            
                            const int idx_cell_wghost_y_B = (i + num_subghosts_velocity[0]) +
                                (j - 1 + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                                (k + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                    subghostcell_dims_velocity[1];
                            
                            const int idx_cell_wghost_y_T = (i + num_subghosts_velocity[0]) +
                                (j + 1 + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                                (k + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                    subghostcell_dims_velocity[1];
                            
                            const int idx_cell_wghost_z_B = (i + num_subghosts_velocity[0]) +
                                (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                                (k - 1 + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                    subghostcell_dims_velocity[1];
                            
                            const int idx_cell_wghost_z_F = (i + num_subghosts_velocity[0]) +
                                (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                                (k + 1 + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                    subghostcell_dims_velocity[1];
                            
                            const int idx_cell_nghost = i +
                                j*interior_dims[0] +
                                k*interior_dims[0]*interior_dims[1];
                            
                            const int idx_face_x_LL = i +
                                (j + 1)*(interior_dims[0] + 3) +
                                (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                            
                            const int idx_face_x_L = (i + 1) +
                                (j + 1)*(interior_dims[0] + 3) +
                                (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                            
                            const int idx_face_x_R = (i + 2) +
                                (j + 1)*(interior_dims[0] + 3) +
                                (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                            
                            const int idx_face_x_RR = (i + 3) +
                                (j + 1)*(interior_dims[0] + 3) +
                                (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                            
                            const int idx_face_y_BB = j +
                                (k + 1)*(interior_dims[1] + 3) +
                                (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                            
                            const int idx_face_y_B = (j + 1) +
                                (k + 1)*(interior_dims[1] + 3) +
                                (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                            
                            const int idx_face_y_T = (j + 2) +
                                (k + 1)*(interior_dims[1] + 3) +
                                (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                            
                            const int idx_face_y_TT = (j + 3) +
                                (k + 1)*(interior_dims[1] + 3) +
                                (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                            
                            const int idx_face_z_BB = k +
                                (i + 1)*(interior_dims[2] + 3) +
                                (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                            
                            const int idx_face_z_B = (k + 1) +
                                (i + 1)*(interior_dims[2] + 3) +
                                (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                            
                            const int idx_face_z_F = (k + 2) +
                                (i + 1)*(interior_dims[2] + 3) +
                                (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                            
                            const int idx_face_z_FF = (k + 3) +
                                (i + 1)*(interior_dims[2] + 3) +
                                (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                            
                            const double& u_LL = velocity_intercell->getPointer(0, 0)[idx_face_x_LL];
                            const double& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                            const double& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                            const double& u_RR = velocity_intercell->getPointer(0, 0)[idx_face_x_RR];
                            
                            const double& v_BB = velocity_intercell->getPointer(1, 1)[idx_face_y_BB];
                            const double& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                            const double& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                            const double& v_TT = velocity_intercell->getPointer(1, 1)[idx_face_y_TT];
                            
                            const double& w_BB = velocity_intercell->getPointer(2, 2)[idx_face_z_BB];
                            const double& w_B = velocity_intercell->getPointer(2, 2)[idx_face_z_B];
                            const double& w_F = velocity_intercell->getPointer(2, 2)[idx_face_z_F];
                            const double& w_FF = velocity_intercell->getPointer(2, 2)[idx_face_z_FF];
                            
                            S[idx_cell_nghost] += dt*Q[ei][idx_cell_wghost]*(
                                (3.0/2*(u_R - u_L) - 3.0/10*(u[idx_cell_wghost_x_R] - u[idx_cell_wghost_x_L]) +
                                 1.0/30*(u_RR - u_LL))/dx[0] +
                                (3.0/2*(v_T - v_B) - 3.0/10*(v[idx_cell_wghost_y_T] - v[idx_cell_wghost_y_B]) +
                                 1.0/30*(v_TT - v_BB))/dx[1] +
                                (3.0/2*(w_F - w_B) - 3.0/10*(w[idx_cell_wghost_z_F] - w[idx_cell_wghost_z_B]) +
                                 1.0/30*(w_FF - w_BB))/dx[2]);
                        }
                    }
                }
            }
        }
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(3))
}




/*
 * Transform physcial variables into characteristic variables.
 * W_array: Characteristic variables.
 * U_array: Physical variables.
 * R_inv_intercell: Projection matrix.
 */
void
ConvectiveFluxReconstructorWCNS56::projectPhysicalVariablesToCharacteristicFields(
    boost::multi_array<double, 2>& W_array,
    const boost::multi_array<const double*, 2>& U_array,
    const boost::multi_array<double, 2>& R_inv_intercell)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(W_array.shape()[0] == U_array.shape()[0]);
    TBOX_ASSERT(W_array.shape()[1] == U_array.shape()[1]);
    TBOX_ASSERT(R_inv_intercell.shape()[0] == R_inv_intercell.shape()[1]);
    TBOX_ASSERT(R_inv_intercell.shape()[0] == U_array.shape()[1]);
#endif
    
    for (int ri = 0; ri < static_cast<int>(R_inv_intercell.shape()[0]); ri++)
    {
        for (int m = 0; m < static_cast<int>(W_array.shape()[0]); m++)
        {
            W_array[m][ri] = 0.0;
            for (int ci = 0; ci < static_cast<int>(R_inv_intercell.shape()[1]); ci++)
            {
                W_array[m][ri] += R_inv_intercell[ri][ci]*(*U_array[m][ci]);
            }
            
        }
    }
}


/*
 * Transform characteristic variables into primitive variables.
 * U: Physcial variables.
 * W: Characteristic variables.
 * R_intercell: Inverse of projection matrix.
 */
void
ConvectiveFluxReconstructorWCNS56::projectCharacteristicVariablesToPhysicalFields(
    std::vector<double>& U,
    const std::vector<double>& W,
    const boost::multi_array<double, 2>& R_intercell)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(U.size() == W.size());
    TBOX_ASSERT(R_intercell.shape()[0] == R_intercell.shape()[1]);
    TBOX_ASSERT(static_cast<int>(R_intercell.shape()[0]) == static_cast<int>(U.size()));
#endif
    
    for (int ri = 0; ri < static_cast<int>(R_intercell.shape()[0]); ri++)
    {
        U[ri] = 0.0;
        for (int ci = 0; ci < static_cast<int>(R_intercell.shape()[1]); ci++)
        {
            U[ri] += R_intercell[ri][ci]*W[ci];
        }
    }
}
