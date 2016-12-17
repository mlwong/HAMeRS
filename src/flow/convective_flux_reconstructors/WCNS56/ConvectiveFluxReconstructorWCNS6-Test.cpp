#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS6-Test.hpp"

#define EPSILON 1e-40

/*
 * Timers interspersed throughout the class.
 */

boost::shared_ptr<tbox::Timer> ConvectiveFluxReconstructorWCNS6_Test::t_characteristic_decomposition;
boost::shared_ptr<tbox::Timer> ConvectiveFluxReconstructorWCNS6_Test::t_WENO_interpolation;
boost::shared_ptr<tbox::Timer> ConvectiveFluxReconstructorWCNS6_Test::t_Riemann_solver;
boost::shared_ptr<tbox::Timer> ConvectiveFluxReconstructorWCNS6_Test::t_reconstruct_flux;

ConvectiveFluxReconstructorWCNS6_Test::ConvectiveFluxReconstructorWCNS6_Test(
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
            convective_flux_reconstructor_db),
        W_array(
            boost::extents[6][d_num_eqn],
            boost::fortran_storage_order()),
        W_minus(d_num_eqn),
        W_plus(d_num_eqn),
        beta(4),
        beta_tilde(4)
{
    d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*4;
    
    d_constant_C = d_convective_flux_reconstructor_db->getDoubleWithDefault("constant_C", 1.0e9);
    d_constant_C = d_convective_flux_reconstructor_db->getDoubleWithDefault("d_constant_C", d_constant_C);
    
    d_constant_p = d_convective_flux_reconstructor_db->getIntegerWithDefault("constant_p", 2);
    d_constant_p = d_convective_flux_reconstructor_db->getIntegerWithDefault("d_constant_p", d_constant_p);
    
    d_constant_q = d_convective_flux_reconstructor_db->getIntegerWithDefault("constant_q", 4);
    d_constant_q = d_convective_flux_reconstructor_db->getIntegerWithDefault("d_constant_q", d_constant_q);
    
    d_constant_alpha_beta = d_convective_flux_reconstructor_db->getDoubleWithDefault("constant_alpha_beta", 35.0);
    d_constant_alpha_beta = d_convective_flux_reconstructor_db->getDoubleWithDefault("d_constant_alpha_beta", d_constant_alpha_beta);
    
    d_weights_c.resize(boost::extents[4][3]);
    d_weights_c[0][0] = 3.0/8;
    d_weights_c[0][1] = -5.0/4;
    d_weights_c[0][2] = 15.0/8;
    d_weights_c[1][0] = -1.0/8;
    d_weights_c[1][1] = 3.0/4;
    d_weights_c[1][2] = 3.0/8;
    d_weights_c[2][0] = 3.0/8;
    d_weights_c[2][1] = 3.0/4;
    d_weights_c[2][2] = -1.0/8;
    d_weights_c[3][0] = 15.0/8;
    d_weights_c[3][1] = -5.0/4;
    d_weights_c[3][2] = 3.0/8;
    
    t_characteristic_decomposition = tbox::TimerManager::getManager()->
        getTimer("ConvectiveFluxReconstructorWCNS6_Test::t_characteristic_decomposition");
    
    t_WENO_interpolation = tbox::TimerManager::getManager()->
        getTimer("ConvectiveFluxReconstructorWCNS6_Test::t_WENO_interpolation");
    
    t_Riemann_solver = tbox::TimerManager::getManager()->
        getTimer("ConvectiveFluxReconstructorWCNS6_Test::t_Riemann_solver");
    
    t_reconstruct_flux = tbox::TimerManager::getManager()->
        getTimer("ConvectiveFluxReconstructorWCNS6_Test::t_reconstruct_flux");
    
}


ConvectiveFluxReconstructorWCNS6_Test::~ConvectiveFluxReconstructorWCNS6_Test()
{
    t_characteristic_decomposition.reset();
    t_WENO_interpolation.reset();
    t_Riemann_solver.reset();
    t_reconstruct_flux.reset();
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorWCNS6_Test object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorWCNS6_Test: this = "
       << (ConvectiveFluxReconstructorWCNS6_Test *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_constant_C = "
       << d_constant_C
       << std::endl;
    os << "d_constant_p = "
       << d_constant_p
       << std::endl;
    os << "d_constant_q = "
       << d_constant_q
       << std::endl;
    os << "d_constant_alpha_beta = "
       << d_constant_alpha_beta
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putDouble("d_constant_C", d_constant_C);
    restart_db->putInteger("d_constant_p", d_constant_p);
    restart_db->putInteger("d_constant_q", d_constant_q);
    restart_db->putDouble("d_constant_alpha_beta", d_constant_alpha_beta);
}


/*
 * Compute the convective fluxes and sources due to hyperbolization
 * of the equations.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::computeConvectiveFluxesAndSources(
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
boost::shared_ptr<pdat::SideData<double> > velocity_intercell(
    new pdat::SideData<double>(interior_box, d_dim.getValue(), hier::IntVector::getOne(d_dim)));
    
    boost::shared_ptr<pdat::CellData<double> > vorticity_magnitude(
        new pdat::CellData<double>(interior_box, 1, d_num_conv_ghosts));
    
boost::shared_ptr<pdat::SideData<double> > convective_flux_midpoint(
    new pdat::SideData<double>(interior_box, d_num_eqn, hier::IntVector::getOne(d_dim)));
    
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
            AVERAGING::SIMPLE);
        
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
                DIRECTION::X_DIRECTION);
            
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
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            
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
            AVERAGING::SIMPLE);
        
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
 * Declare temporary data containers.
 */

std::vector<boost::shared_ptr<pdat::SideData<double> > > projection_variables;
std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > > characteristic_variables;

std::vector<boost::shared_ptr<pdat::SideData<double> > > characteristic_variables_minus;
std::vector<boost::shared_ptr<pdat::SideData<double> > > characteristic_variables_plus;

std::vector<boost::shared_ptr<pdat::SideData<double> > > primitive_variables_minus;
std::vector<boost::shared_ptr<pdat::SideData<double> > > primitive_variables_plus;

/*
 * Initialize temporary data containers.
 */

int num_projection_var = d_flow_model->getNumberOfProjectionVariablesForPrimitiveVariables();
projection_variables.reserve(num_projection_var);
for (int vi = 0; vi < num_projection_var; vi++)
{
    projection_variables.push_back(boost::make_shared<pdat::SideData<double> >(
        interior_box, 1, hier::IntVector::getOne(d_dim)));
}

characteristic_variables.resize(6);
for (int m = 0; m < 6; m++)
{
    characteristic_variables[m].reserve(d_num_eqn);
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        characteristic_variables[m].push_back(boost::make_shared<pdat::SideData<double> >(
            interior_box, 1, hier::IntVector::getOne(d_dim)));
    }
}

characteristic_variables_minus.reserve(d_num_eqn);
characteristic_variables_plus.reserve(d_num_eqn);
primitive_variables_minus.reserve(d_num_eqn);
primitive_variables_plus.reserve(d_num_eqn);
for (int ei = 0; ei < d_num_eqn; ei++)
{
    characteristic_variables_minus.push_back(boost::make_shared<pdat::SideData<double> >(
        interior_box, 1, hier::IntVector::getOne(d_dim)));
    
    characteristic_variables_plus.push_back(boost::make_shared<pdat::SideData<double> >(
        interior_box, 1, hier::IntVector::getOne(d_dim)));
    
    primitive_variables_minus.push_back(boost::make_shared<pdat::SideData<double> >(
        interior_box, 1, hier::IntVector::getOne(d_dim)));
    
    primitive_variables_plus.push_back(boost::make_shared<pdat::SideData<double> >(
        interior_box, 1, hier::IntVector::getOne(d_dim)));
}

/*
 * Compute global side data of the projection variables for transformation between
 * primitive variables and characteristic variables.
 */
t_characteristic_decomposition->start();

d_flow_model->computeGlobalSideDataProjectionVariablesForPrimitiveVariables(
    projection_variables);

t_characteristic_decomposition->stop();

/*
 * Transform primitive variables to characteristic variables.
 */

t_characteristic_decomposition->start();

for (int m = 0; m < 6; m++)
{
    d_flow_model->computeGlobalSideDataCharacteristicVariablesFromPrimitiveVariables(
        characteristic_variables[m],
        primitive_variables,
        projection_variables,
        m - 3);
}

t_characteristic_decomposition->stop();

/*
 * Peform WENO interpolation in the x-direction.
 */

t_WENO_interpolation->start();

std::vector<double*> W_array_x;
W_array_x.resize(6);

for (int ei = 0; ei < d_num_eqn; ei++)
{
    for (int m = 0; m < 6; m++)
    {
        W_array_x[m] = characteristic_variables[m][ei]->getPointer(0);
    }
    
    double* W_L_ptr = characteristic_variables_minus[ei]->getPointer(0);
    double* W_R_ptr = characteristic_variables_plus[ei]->getPointer(0);
    
    const int interior_dim_0 = interior_dims[0];
    const int interior_dim_1 = interior_dims[1];
    
    #pragma ivdep
    for (int j = 0; j < interior_dim_1; j++)
    {
        for (int i = -1; i < interior_dim_0 + 2; i++)
        {
            // Compute the linear index of the side.
            const int idx_side_x = (i + 1) +
                (j + 1)*(interior_dim_0 + 3);
            
            #pragma forceinline
            performWENOInterpolation_new(
                W_L_ptr,
                W_R_ptr,
                W_array_x,
                idx_side_x);
        }
    }
}

t_WENO_interpolation->stop();

/*
 * Peform WENO interpolation in the y-direction.
 */

t_WENO_interpolation->start();

std::vector<double*> W_array_y;
W_array_y.resize(6);

for (int ei = 0; ei < d_num_eqn; ei++)
{
    for (int m = 0; m < 6; m++)
    {
        W_array_y[m] = characteristic_variables[m][ei]->getPointer(1);
    }
    
    double* W_B_ptr = characteristic_variables_minus[ei]->getPointer(1);
    double* W_T_ptr = characteristic_variables_plus[ei]->getPointer(1);
    
    const int interior_dim_0 = interior_dims[0];
    const int interior_dim_1 = interior_dims[1];
    
    #pragma ivdep
    for (int j = -1; j < interior_dim_1 + 2; j++)
    {
        for (int i = 0; i < interior_dim_0; i++)
        {
            // Compute the linear index of the side.
            const int idx_side_y = (i + 1) +
                (j + 1)*(interior_dim_0 + 2);
            
            #pragma forceinline
            performWENOInterpolation_new(
                W_B_ptr,
                W_T_ptr,
                W_array_y,
                idx_side_y);
        }
    }
}

t_WENO_interpolation->stop();

/*
 * Transform characteristic variables back to primitive variables.
 */

t_characteristic_decomposition->start();

d_flow_model->computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables(
    primitive_variables_minus,
    characteristic_variables_minus,
    projection_variables);

d_flow_model->computeGlobalSideDataPrimitiveVariablesFromCharacteristicVariables(
    primitive_variables_plus,
    characteristic_variables_plus,
    projection_variables);

t_characteristic_decomposition->stop();

/*
 * Compute mid-point flux in the x-direction.
 */

t_Riemann_solver->start();

std::vector<double*> V_x_L;
std::vector<double*> V_x_R;
V_x_L.resize(d_num_eqn);
V_x_R.resize(d_num_eqn);

for (int ei = 0; ei < d_num_eqn; ei++)
{
    V_x_L[ei] = primitive_variables_minus[ei]->getPointer(0);
    V_x_R[ei] = primitive_variables_plus[ei]->getPointer(0);
}

// Declare and initialize containers to store the references to the
// WENO interpolated values.
std::vector<boost::reference_wrapper<double> > V_x_L_ref;
std::vector<boost::reference_wrapper<double> > V_x_R_ref;
V_x_L_ref.reserve(d_num_eqn);
V_x_R_ref.reserve(d_num_eqn);

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
        // Compute the linear index of the side.
        const int idx_side_x = (i + 1) +
            (j + 1)*(interior_dims[0] + 3);
            
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
            F_x_midpoint_ref.push_back(boost::ref(F_x_midpoint[ei][idx_side_x]));
        }
        
        // Initialize container that stores the references to the mid-point velocity.
        for (int di = 0; di < d_dim.getValue(); di++)
        {
            vel_x_midpoint_ref.push_back(
                boost::ref(velocity_intercell->getPointer(0, di)[idx_side_x]));
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_x_L_ref.push_back(boost::ref(V_x_L[ei][idx_side_x]));
            V_x_R_ref.push_back(boost::ref(V_x_R[ei][idx_side_x]));
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
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC_HLL);
        }
        else
        {
            d_flow_model->
                computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                    F_x_midpoint_ref,
                    vel_x_midpoint_ref,
                    V_x_L_ref,
                    V_x_R_ref,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
        }
        
        F_x_midpoint_ref.clear();
        vel_x_midpoint_ref.clear();
        V_x_L_ref.clear();
        V_x_R_ref.clear();
    }
}

t_Riemann_solver->stop();

/*
 * Compute mid-point flux in the y-direction.
 */

t_Riemann_solver->start();

std::vector<double*> V_y_B;
std::vector<double*> V_y_T;
V_y_B.resize(d_num_eqn);
V_y_T.resize(d_num_eqn);

for (int ei = 0; ei < d_num_eqn; ei++)
{
    V_y_B[ei] = primitive_variables_minus[ei]->getPointer(1);
    V_y_T[ei] = primitive_variables_plus[ei]->getPointer(1);
}

// Declare and initialize containers to store the references to the
// WENO interpolated values.
std::vector<boost::reference_wrapper<double> > V_y_B_ref;
std::vector<boost::reference_wrapper<double> > V_y_T_ref;
V_y_B_ref.reserve(d_num_eqn);
V_y_T_ref.reserve(d_num_eqn);

// Declare container to store the references to the mid-point flux.
std::vector<boost::reference_wrapper<double> > F_y_midpoint_ref;
F_y_midpoint_ref.reserve(d_num_eqn);

// Declare container to store the references to the mid-point velocity.
std::vector<boost::reference_wrapper<double> > vel_y_midpoint_ref;
vel_y_midpoint_ref.reserve(d_dim.getValue());

for (int j = -1; j < interior_dims[1] + 2; j++)
{
    for (int i = 0; i < interior_dims[0]; i++)
    {
        // Compute the linear index of the side.
        const int idx_side_y = (i + 1) +
            (j + 1)*(interior_dims[0] + 2);
        
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
            F_y_midpoint_ref.push_back(boost::ref(F_y_midpoint[ei][idx_side_y]));
        }
        
        // Initialize container that stores the references to the mid-point velocity.
        for (int di = 0; di < d_dim.getValue(); di++)
        {
            vel_y_midpoint_ref.push_back(
                boost::ref(velocity_intercell->getPointer(1, di)[idx_side_y]));
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_y_B_ref.push_back(boost::ref(V_y_B[ei][idx_side_y]));
            V_y_T_ref.push_back(boost::ref(V_y_T[ei][idx_side_y]));
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
                    DIRECTION::Y_DIRECTION,
                    RIEMANN_SOLVER::HLLC_HLL);
        }
        else
        {
            d_flow_model->
                computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                    F_y_midpoint_ref,
                    vel_y_midpoint_ref,
                    V_y_B_ref,
                    V_y_T_ref,
                    DIRECTION::Y_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
        }
        
        F_y_midpoint_ref.clear();
        vel_y_midpoint_ref.clear();
        V_y_B_ref.clear();
        V_y_T_ref.clear();
    }
}

t_Riemann_solver->stop();

        /*
         * Reconstruct the flux in the x direction.
         */
        
t_reconstruct_flux->start();
        
        for (int j = 0; j < interior_dims[1]; j++)
        {
            for (int i = 0; i < interior_dims[0] + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i +
                    j*(interior_dims[0] + 1);
                
const int idx_midpoint_x = (i + 1) +
    (j + 1)*(interior_dims[0] + 3);

const int idx_midpoint_x_L = i +
    (j + 1)*(interior_dims[0] + 3);

const int idx_midpoint_x_R = (i + 2) +
    (j + 1)*(interior_dims[0] + 3);
                
                const int idx_node_L = (i - 1 + num_subghosts_convective_flux_x[0]) +
                    (j + num_subghosts_convective_flux_x[1])*subghostcell_dims_convective_flux_x[0];
                
                const int idx_node_R = (i + num_subghosts_convective_flux_x[0]) +
                    (j + num_subghosts_convective_flux_x[1])*subghostcell_dims_convective_flux_x[0];
                
                // Compute the fluxes.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    convective_flux->getPointer(0, ei)[idx_face_x] =
                        dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x_R] +
                        F_x_midpoint[ei][idx_midpoint_x_L]) -
                        3.0/10*(F_x_node[ei][idx_node_R] +
                        F_x_node[ei][idx_node_L]) +
                        23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
                }
            }
        }
        
t_reconstruct_flux->stop();
        
        /*
         * Reconstruct the flux in the y direction.
         */
        
t_reconstruct_flux->start();
        
        for (int i = 0; i < interior_dims[0]; i++)
        {
            for (int j = 0; j < interior_dims[1] + 1; j++)
            {
                // Compute the linear indices.
                const int idx_face_y = j +
                    i*(interior_dims[1] + 1);
                
/*
const int idx_midpoint_y = (j + 1) +
    (i + 1)*(interior_dims[1] + 3);
*/
const int idx_midpoint_y = (i + 1) +
    (j + 1)*(interior_dims[0] + 2);

const int idx_midpoint_y_B = (i + 1) +
    j*(interior_dims[0] + 2);

const int idx_midpoint_y_T = (i + 1) +
    (j + 2)*(interior_dims[0] + 2);
                
                const int idx_node_B = (i + num_subghosts_convective_flux_y[0]) +
                    (j - 1 + num_subghosts_convective_flux_y[1])*subghostcell_dims_convective_flux_y[0];
                
                const int idx_node_T = (i + num_subghosts_convective_flux_y[0]) +
                    (j + num_subghosts_convective_flux_y[1])*subghostcell_dims_convective_flux_y[0];
                
                // Compute the fluxes.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    convective_flux->getPointer(1, ei)[idx_face_y] =
                        dt*(1.0/30*(F_y_midpoint[ei][idx_midpoint_y_T] +
                        F_y_midpoint[ei][idx_midpoint_y_B]) -
                        3.0/10*(F_y_node[ei][idx_node_T] +
                        F_y_node[ei][idx_node_B]) +
                        23.0/15*F_y_midpoint[ei][idx_midpoint_y]);
                }
            }
        }
        
t_reconstruct_flux->stop();
        
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
                        
                        const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_velocity[0]) +
                            (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                        
                        const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_velocity[0]) +
                            (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                        
                        const int idx_cell_wghost_y_B = (i + num_subghosts_velocity[0]) +
                            (j - 1 + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                        
                        const int idx_cell_wghost_y_T = (i + num_subghosts_velocity[0]) +
                            (j + 1 + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                        
                        const int idx_cell_nghost = i + j*interior_dims[0];
                        
                        const int idx_side_x_LL = i + 
                            (j + 1)*(interior_dims[0] + 3);
                        
                        const int idx_side_x_L = (i + 1) +
                            (j + 1)*(interior_dims[0] + 3);
                        
                        const int idx_side_x_R = (i + 2) +
                            (j + 1)*(interior_dims[0] + 3);
                        
                        const int idx_side_x_RR = (i + 3) +
                            (j + 1)*(interior_dims[0] + 3);
                        
/*
const int idx_face_y_BB = j +
    (i + 1)*(interior_dims[1] + 3);

const int idx_face_y_B = (j + 1) +
    (i + 1)*(interior_dims[1] + 3);

const int idx_face_y_T = (j + 2) +
    (i + 1)*(interior_dims[1] + 3);

const int idx_face_y_TT = (j + 3) +
    (i + 1)*(interior_dims[1] + 3);
*/

const int idx_side_y_BB = (i + 1) +
    j*(interior_dims[0] + 2);

const int idx_side_y_B = (i + 1) +
    (j + 1)*(interior_dims[0] + 2);

const int idx_side_y_T = (i + 1) +
    (j + 2)*(interior_dims[0] + 2);

const int idx_side_y_TT = (i + 1) +
    (j + 3)*(interior_dims[0] + 2);
                        
                        const double& u_LL = velocity_intercell->getPointer(0, 0)[idx_side_x_LL];
                        const double& u_L = velocity_intercell->getPointer(0, 0)[idx_side_x_L];
                        const double& u_R = velocity_intercell->getPointer(0, 0)[idx_side_x_R];
                        const double& u_RR = velocity_intercell->getPointer(0, 0)[idx_side_x_RR];
                        
                        const double& v_BB = velocity_intercell->getPointer(1, 1)[idx_side_y_BB];
                        const double& v_B = velocity_intercell->getPointer(1, 1)[idx_side_y_B];
                        const double& v_T = velocity_intercell->getPointer(1, 1)[idx_side_y_T];
                        const double& v_TT = velocity_intercell->getPointer(1, 1)[idx_side_y_TT];
                        
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
            AVERAGING::SIMPLE);
        
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
                        DIRECTION::X_DIRECTION);
                    
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
                                DIRECTION::X_DIRECTION,
                                RIEMANN_SOLVER::HLLC_HLL);
                    }
                    else
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_x_midpoint_ref,
                                vel_x_midpoint_ref,
                                V_x_L_ref,
                                V_x_R_ref,
                                DIRECTION::X_DIRECTION,
                                RIEMANN_SOLVER::HLLC);
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
                        DIRECTION::Y_DIRECTION);
                    
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
                                DIRECTION::Y_DIRECTION,
                                RIEMANN_SOLVER::HLLC_HLL);
                    }
                    else
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_y_midpoint_ref,
                                vel_y_midpoint_ref,
                                V_y_B_ref,
                                V_y_T_ref,
                                DIRECTION::Y_DIRECTION,
                                RIEMANN_SOLVER::HLLC);
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
                        DIRECTION::Z_DIRECTION);
                    
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
                                DIRECTION::Z_DIRECTION,
                                RIEMANN_SOLVER::HLLC_HLL);
                    }
                    else
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_z_midpoint_ref,
                                vel_z_midpoint_ref,
                                V_z_B_ref,
                                V_z_F_ref,
                                DIRECTION::Z_DIRECTION,
                                RIEMANN_SOLVER::HLLC);
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
ConvectiveFluxReconstructorWCNS6_Test::projectPhysicalVariablesToCharacteristicFields(
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
ConvectiveFluxReconstructorWCNS6_Test::projectCharacteristicVariablesToPhysicalFields(
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


/*
 * Compute sigma's.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::computeSigma(
    double& sigma,
    const boost::multi_array_ref<double, 2>::const_array_view<1>::type& W_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
#endif
    
    /*
     * Compute the sigma.
     */
    
    const double alpha_1 = W_array[2] - W_array[1];
    const double alpha_2 = W_array[3] - W_array[2];
    const double alpha_3 = W_array[4] - W_array[3];
    
    const double theta_1 = fabs(alpha_1 - alpha_2)/(fabs(alpha_1) + fabs(alpha_2) + EPSILON);
    const double theta_2 = fabs(alpha_2 - alpha_3)/(fabs(alpha_2) + fabs(alpha_3) + EPSILON);
    
    sigma = fmax(theta_1, theta_2);
}


/*
 * Compute beta's.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::computeBeta(
    std::vector<double>& beta,
    const boost::multi_array_ref<double, 2>::const_array_view<1>::type& W_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(beta.size()) == 4);
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
#endif
    
    beta[0] = 1.0/3*(W_array[0]*(4*W_array[0] - 19*W_array[1] + 11*W_array[2]) +
        W_array[1]*(25*W_array[1] - 31*W_array[2]) + 10*W_array[2]*W_array[2]);
    
    beta[1] = 1.0/3*(W_array[1]*(4*W_array[1] - 13*W_array[2] + 5*W_array[3]) +
        13*W_array[2]*(W_array[2] - W_array[3]) + 4*W_array[3]*W_array[3]);
    
    beta[2] = 1.0/3*(W_array[2]*(10*W_array[2] - 31*W_array[3] + 11*W_array[4]) +
        W_array[3]*(25*W_array[3] - 19*W_array[4]) + 4*W_array[4]*W_array[4]);
    
    beta[3] = 1.0/232243200*(W_array[0]*(525910327*W_array[0] - 4562164630*W_array[1] +
        7799501420*W_array[2] - 6610694540*W_array[3] + 2794296070*W_array[4] -
        472758974*W_array[5]) + 5*W_array[1]*(2146987907*W_array[1] - 7722406988*W_array[2] +
        6763559276*W_array[3] - 2926461814*W_array[4] + 503766638*W_array[5]) +
        20*W_array[2]*(1833221603*W_array[2] - 3358664662*W_array[3] + 1495974539*W_array[4] -
        263126407*W_array[5]) + 20*W_array[3]*(1607794163*W_array[3] - 1486026707*W_array[4] +
        268747951*W_array[5]) +  5*W_array[4]*(1432381427*W_array[4] - 536951582*W_array[5]) +
        263126407*W_array[5]*W_array[5]);
}


/*
 * Compute beta_tilde's.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::computeBetaTilde(
    std::vector<double>& beta_tilde,
    const boost::multi_array_ref<double, 2>::const_array_view<1>::type& W_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(beta_tilde.size()) == 4);
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
#endif
    
    beta_tilde[0] = 1.0/3*(W_array[5]*(4*W_array[5] - 19*W_array[4] + 11*W_array[3]) +
        W_array[4]*(25*W_array[4] - 31*W_array[3]) + 10*W_array[3]*W_array[3]);
    
    beta_tilde[1] = 1.0/3*(W_array[4]*(4*W_array[4] - 13*W_array[3] + 5*W_array[2]) +
        13*W_array[3]*(W_array[3] - W_array[2]) + 4*W_array[2]*W_array[2]);
    
    beta_tilde[2] = 1.0/3*(W_array[3]*(10*W_array[3] - 31*W_array[2] + 11*W_array[1]) +
        W_array[2]*(25*W_array[2] - 19*W_array[1]) + 4*W_array[1]*W_array[1]);
    
    beta_tilde[3] = 1.0/232243200*(W_array[5]*(525910327*W_array[5] - 4562164630*W_array[4] +
        7799501420*W_array[3] - 6610694540*W_array[2] + 2794296070*W_array[1] -
        472758974*W_array[0]) + 5*W_array[4]*(2146987907*W_array[4] - 7722406988*W_array[3] +
        6763559276*W_array[2] - 2926461814*W_array[1] + 503766638*W_array[0]) +
        20*W_array[3]*(1833221603*W_array[3] - 3358664662*W_array[2] + 1495974539*W_array[1] -
        263126407*W_array[0]) + 20*W_array[2]*(1607794163*W_array[2] -
        1486026707*W_array[1] + 268747951*W_array[0]) + 5*W_array[1]*(1432381427*W_array[1] -
        536951582*W_array[0])+263126407*W_array[0]*W_array[0]);
}


/*
 * Compute sigma's.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::computeSigma(
    double& sigma,
    const boost::multi_array_ref<const double*, 2>::const_array_view<1>::type& U_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
#endif
    
    /*
     * Compute the sigma.
     */
    
    const double alpha_1 = *(U_array[2]) - *(U_array[1]);
    const double alpha_2 = *(U_array[3]) - *(U_array[2]);
    const double alpha_3 = *(U_array[4]) - *(U_array[3]);
    
    const double theta_1 = fabs(alpha_1 - alpha_2)/(fabs(alpha_1) + fabs(alpha_2) + EPSILON);
    const double theta_2 = fabs(alpha_2 - alpha_3)/(fabs(alpha_2) + fabs(alpha_3) + EPSILON);
    
    sigma = fmax(theta_1, theta_2);
}


/*
 * Compute beta's.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::computeBeta(
    std::vector<double>& beta,
    const boost::multi_array_ref<const double*, 2>::const_array_view<1>::type& U_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(beta.size()) == 4);
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
#endif
    
    beta[0] = 1.0/3*((*U_array[0])*(4*(*U_array[0]) - 19*(*U_array[1]) + 11*(*U_array[2])) +
        (*U_array[1])*(25*(*U_array[1]) - 31*(*U_array[2])) + 10*(*U_array[2])*(*U_array[2]));
    
    beta[1] = 1.0/3*((*U_array[1])*(4*(*U_array[1]) - 13*(*U_array[2]) + 5*(*U_array[3])) +
        13*(*U_array[2])*((*U_array[2]) - (*U_array[3])) + 4*(*U_array[3])*(*U_array[3]));
    
    beta[2] = 1.0/3*((*U_array[2])*(10*(*U_array[2]) - 31*(*U_array[3]) + 11*(*U_array[4])) +
        (*U_array[3])*(25*(*U_array[3]) - 19*(*U_array[4])) + 4*(*U_array[4])*(*U_array[4]));
    
    beta[3] = 1.0/232243200*((*U_array[0])*(525910327*(*U_array[0]) - 4562164630*(*U_array[1]) +
        7799501420*(*U_array[2]) - 6610694540*(*U_array[3]) + 2794296070*(*U_array[4]) -
        472758974*(*U_array[5])) + 5*(*U_array[1])*(2146987907*(*U_array[1]) - 7722406988*(*U_array[2]) +
        6763559276*(*U_array[3]) - 2926461814*(*U_array[4]) + 503766638*(*U_array[5])) +
        20*(*U_array[2])*(1833221603*(*U_array[2]) - 3358664662*(*U_array[3]) + 1495974539*(*U_array[4]) -
        263126407*(*U_array[5])) + 20*(*U_array[3])*(1607794163*(*U_array[3]) - 1486026707*(*U_array[4]) +
        268747951*(*U_array[5])) +  5*(*U_array[4])*(1432381427*(*U_array[4]) - 536951582*(*U_array[5])) +
        263126407*(*U_array[5])*(*U_array[5]));
}


/*
 * Compute beta_tilde's.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::computeBetaTilde(
    std::vector<double>& beta_tilde,
    const boost::multi_array_ref<const double*, 2>::const_array_view<1>::type& U_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(beta_tilde.size()) == 4);
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
#endif
    
    beta_tilde[0] = 1.0/3*((*U_array[5])*(4*(*U_array[5]) - 19*(*U_array[4]) + 11*(*U_array[3])) +
        (*U_array[4])*(25*(*U_array[4]) - 31*(*U_array[3])) + 10*(*U_array[3])*(*U_array[3]));
    
    beta_tilde[1] = 1.0/3*((*U_array[4])*(4*(*U_array[4]) - 13*(*U_array[3]) + 5*(*U_array[2])) +
        13*(*U_array[3])*((*U_array[3]) - (*U_array[2])) + 4*(*U_array[2])*(*U_array[2]));
    
    beta_tilde[2] = 1.0/3*((*U_array[3])*(10*(*U_array[3]) - 31*(*U_array[2]) + 11*(*U_array[1])) +
        (*U_array[2])*(25*(*U_array[2]) - 19*(*U_array[1])) + 4*(*U_array[1])*(*U_array[1]));
    
    beta_tilde[3] = 1.0/232243200*((*U_array[5])*(525910327*(*U_array[5]) - 4562164630*(*U_array[4]) +
        7799501420*(*U_array[3]) - 6610694540*(*U_array[2]) + 2794296070*(*U_array[1]) -
        472758974*(*U_array[0])) + 5*(*U_array[4])*(2146987907*(*U_array[4]) - 7722406988*(*U_array[3]) +
        6763559276*(*U_array[2]) - 2926461814*(*U_array[1]) + 503766638*(*U_array[0])) +
        20*(*U_array[3])*(1833221603*(*U_array[3]) - 3358664662*(*U_array[2]) + 1495974539*(*U_array[1]) -
        263126407*(*U_array[0])) + 20*(*U_array[2])*(1607794163*(*U_array[2]) -
        1486026707*(*U_array[1]) + 268747951*(*U_array[0])) + 5*(*U_array[1])*(1432381427*(*U_array[1]) -
        536951582*(*U_array[0]))+263126407*(*U_array[0])*(*U_array[0]));
}


/*
 * Compute sigma's.
 */
inline void
ConvectiveFluxReconstructorWCNS6_Test::computeSigma(
    double& sigma,
    const std::vector<double*>& U_array,
    const int& idx_face)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.size()) == 6);
#endif
    
    /*
     * Compute the sigma.
     */
    
    const double alpha_1 = U_array[2][idx_face] - U_array[1][idx_face];
    const double alpha_2 = U_array[3][idx_face] - U_array[2][idx_face];
    const double alpha_3 = U_array[4][idx_face] - U_array[3][idx_face];
    
    const double theta_1 = fabs(alpha_1 - alpha_2)/(fabs(alpha_1) + fabs(alpha_2) + EPSILON);
    const double theta_2 = fabs(alpha_2 - alpha_3)/(fabs(alpha_2) + fabs(alpha_3) + EPSILON);
    
    sigma = fmax(theta_1, theta_2);
}


/*
 * Compute beta's.
 */
inline void
ConvectiveFluxReconstructorWCNS6_Test::computeBeta(
    std::vector<double>& beta,
    const std::vector<double*>& U_array,
    const int& idx_face)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(beta.size()) == 4);
    TBOX_ASSERT(static_cast<int>(W_array.size()) == 6);
#endif
    
    beta[0] = 1.0/3*(U_array[0][idx_face]*(4*U_array[0][idx_face] - 19*U_array[1][idx_face] + 11*U_array[2][idx_face]) +
        U_array[1][idx_face]*(25*U_array[1][idx_face] - 31*U_array[2][idx_face]) + 10*U_array[2][idx_face]*U_array[2][idx_face]);
    
    beta[1] = 1.0/3*(U_array[1][idx_face]*(4*U_array[1][idx_face] - 13*U_array[2][idx_face] + 5*U_array[3][idx_face]) +
        13*U_array[2][idx_face]*(U_array[2][idx_face] - U_array[3][idx_face]) + 4*U_array[3][idx_face]*U_array[3][idx_face]);
    
    beta[2] = 1.0/3*(U_array[2][idx_face]*(10*U_array[2][idx_face] - 31*U_array[3][idx_face] + 11*U_array[4][idx_face]) +
        U_array[3][idx_face]*(25*U_array[3][idx_face] - 19*U_array[4][idx_face]) + 4*U_array[4][idx_face]*U_array[4][idx_face]);
    
    beta[3] = 1.0/232243200*(U_array[0][idx_face]*(525910327*U_array[0][idx_face] - 4562164630*U_array[1][idx_face] +
        7799501420*U_array[2][idx_face] - 6610694540*U_array[3][idx_face] + 2794296070*U_array[4][idx_face] -
        472758974*U_array[5][idx_face]) + 5*U_array[1][idx_face]*(2146987907*U_array[1][idx_face] - 7722406988*U_array[2][idx_face] +
        6763559276*U_array[3][idx_face] - 2926461814*U_array[4][idx_face] + 503766638*U_array[5][idx_face]) +
        20*U_array[2][idx_face]*(1833221603*U_array[2][idx_face] - 3358664662*U_array[3][idx_face] + 1495974539*U_array[4][idx_face] -
        263126407*U_array[5][idx_face]) + 20*U_array[3][idx_face]*(1607794163*U_array[3][idx_face] - 1486026707*U_array[4][idx_face] +
        268747951*U_array[5][idx_face]) +  5*U_array[4][idx_face]*(1432381427*U_array[4][idx_face] - 536951582*U_array[5][idx_face]) +
        263126407*U_array[5][idx_face]*U_array[5][idx_face]);
}


/*
 * Compute beta_tilde's.
 */
inline void
ConvectiveFluxReconstructorWCNS6_Test::computeBetaTilde(
    std::vector<double>& beta_tilde,
    const std::vector<double*>& U_array,
    const int& idx_face)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(beta_tilde.size()) == 4);
    TBOX_ASSERT(static_cast<int>(W_array.size()) == 6);
#endif
    
    beta_tilde[0] = 1.0/3*(U_array[5][idx_face]*(4*U_array[5][idx_face] - 19*U_array[4][idx_face] + 11*U_array[3][idx_face]) +
        U_array[4][idx_face]*(25*U_array[4][idx_face] - 31*U_array[3][idx_face]) + 10*U_array[3][idx_face]*U_array[3][idx_face]);
    
    beta_tilde[1] = 1.0/3*(U_array[4][idx_face]*(4*U_array[4][idx_face] - 13*U_array[3][idx_face] + 5*U_array[2][idx_face]) +
        13*U_array[3][idx_face]*(U_array[3][idx_face] - U_array[2][idx_face]) + 4*U_array[2][idx_face]*U_array[2][idx_face]);
    
    beta_tilde[2] = 1.0/3*(U_array[3][idx_face]*(10*U_array[3][idx_face] - 31*U_array[2][idx_face] + 11*U_array[1][idx_face]) +
        U_array[2][idx_face]*(25*U_array[2][idx_face] - 19*U_array[1][idx_face]) + 4*U_array[1][idx_face]*U_array[1][idx_face]);
    
    beta_tilde[3] = 1.0/232243200*(U_array[5][idx_face]*(525910327*U_array[5][idx_face] - 4562164630*U_array[4][idx_face] +
        7799501420*U_array[3][idx_face] - 6610694540*U_array[2][idx_face] + 2794296070*U_array[1][idx_face] -
        472758974*U_array[0][idx_face]) + 5*U_array[4][idx_face]*(2146987907*U_array[4][idx_face] - 7722406988*U_array[3][idx_face] +
        6763559276*U_array[2][idx_face] - 2926461814*U_array[1][idx_face] + 503766638*U_array[0][idx_face]) +
        20*U_array[3][idx_face]*(1833221603*U_array[3][idx_face] - 3358664662*U_array[2][idx_face] + 1495974539*U_array[1][idx_face] -
        263126407*U_array[0][idx_face]) + 20*U_array[2][idx_face]*(1607794163*U_array[2][idx_face] -
        1486026707*U_array[1][idx_face] + 268747951*U_array[0][idx_face]) + 5*U_array[1][idx_face]*(1432381427*U_array[1][idx_face] -
        536951582*U_array[0][idx_face])+263126407*U_array[0][idx_face]*U_array[0][idx_face]);
}


/*
 * Compute sigma's.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::computeSigma(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_sigma,
    const std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > >& variables_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(variables_sigma.size()) == d_num_eqn);
    TBOX_ASSERT(variables_sigma[0]->getGhostCellWidth == hier::IntVector::getOne(d_dim));
#endif
    
    /*
     * Get the dimensions of interior box, box including ghost cells and number of ghost cells.
     */
    
    const hier::IntVector interior_dims = variables_sigma[0]->getBox().numberCells();
    const hier::IntVector ghostcell_dims = variables_sigma[0]->getGhostBox().numberCells();
    const hier::IntVector num_ghosts = variables_sigma[0]->getGhostCellWidth();
    
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ERROR(variables_sigma[ei]->getBox().numberCells() == interior_dims);
        TBOX_ERROR(variables_sigma[ei]->getGhostBox().numberCells() == ghostcell_dims);
    }
    
    TBOX_ASSERT(static_cast<int>(variables_array.size()) == 6);
    
    for (int m = 0; m < 6; m++)
    {
        TBOX_ASSERT(static_cast<int>(variables_array[m].size()) == d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ERROR(variables_array[m][ei]->getBox().numberCells() == interior_dims);
            TBOX_ERROR(variables_array[m][ei]->getGhostBox().numberCells() == ghostcell_dims);
        }
    }
#endif
    
    /*
     * Declare the pointers to the sigma's.
     */
    
    std::vector<double*> sigma;
    sigma.resize(d_num_eqn);
    
    /*
     * Declare the pointers to the array of variables.
     */
    
    std::vector<double*> U_array;
    U_array.resize(6);
    
    /*
     * Compute the sigma's.
     */
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        const int num_ghosts_0 = num_ghosts[0];
        
        /*
         * Compute sigma in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            sigma[ei] = variables_sigma[ei]->getPointer(0);
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(0);
            }
            
            for (int i = -num_ghosts_0;
                 i < interior_dim_0 + 1 + num_ghosts_0;
                 i++)
            {
                // Compute the linear index.
                const int idx_side = i + num_ghosts_0;
                
                double alpha_1 = U_array[2][idx_side] - U_array[1][idx_side];
                double alpha_2 = U_array[3][idx_side] - U_array[2][idx_side];
                double alpha_3 = U_array[4][idx_side] - U_array[3][idx_side];
                
                double theta_1 = fabs(alpha_1 - alpha_2)/(fabs(alpha_1) + fabs(alpha_2) + EPSILON);
                double theta_2 = fabs(alpha_2 - alpha_3)/(fabs(alpha_2) + fabs(alpha_3) + EPSILON);
                
                sigma[ei][idx_side] = fmax(theta_1, theta_2);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int ghostcell_dim_0 = ghostcell_dims[0];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        
        /*
         * Compute sigma in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            sigma[ei] = variables_sigma[ei]->getPointer(0);
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(0);
            }
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                for (int i = -num_ghosts_0;
                     i < interior_dim_1 + 1 + num_ghosts_0;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_side = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*(ghostcell_dim_0 + 1);
                    
                    double alpha_1 = U_array[2][idx_side] - U_array[1][idx_side];
                    double alpha_2 = U_array[3][idx_side] - U_array[2][idx_side];
                    double alpha_3 = U_array[4][idx_side] - U_array[3][idx_side];
                    
                    double theta_1 = fabs(alpha_1 - alpha_2)/(fabs(alpha_1) + fabs(alpha_2) + EPSILON);
                    double theta_2 = fabs(alpha_2 - alpha_3)/(fabs(alpha_2) + fabs(alpha_3) + EPSILON);
                    
                    sigma[ei][idx_side] = fmax(theta_1, theta_2);
                }
            }
        }
        
        /*
         * Compute sigma in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            sigma[ei] = variables_sigma[ei]->getPointer(1);
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(1);
            }
            
            for (int j = -num_ghosts_1;
                 j < interior_dim_1 + 1 + num_ghosts_1;
                 j++)
            {
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_side = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*ghostcell_dim_0;
                    
                    double alpha_1 = U_array[2][idx_side] - U_array[1][idx_side];
                    double alpha_2 = U_array[3][idx_side] - U_array[2][idx_side];
                    double alpha_3 = U_array[4][idx_side] - U_array[3][idx_side];
                    
                    double theta_1 = fabs(alpha_1 - alpha_2)/(fabs(alpha_1) + fabs(alpha_2) + EPSILON);
                    double theta_2 = fabs(alpha_2 - alpha_3)/(fabs(alpha_2) + fabs(alpha_3) + EPSILON);
                    
                    sigma[ei][idx_side] = fmax(theta_1, theta_2);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int ghostcell_dim_0 = ghostcell_dims[0];
        const int ghostcell_dim_1 = ghostcell_dims[1];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        const int num_ghosts_2 = num_ghosts[2];
        
        /*
         * Compute sigma in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            sigma[ei] = variables_sigma[ei]->getPointer(0);
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(0);
            }
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    for (int i = -num_ghosts_0;
                         i < interior_dim_1 + 1 + num_ghosts_0;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx_side = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*(ghostcell_dim_0 + 1) +
                            (k + num_ghosts_2)*(ghostcell_dim_0 + 1)*
                                ghostcell_dim_1;
                        
                        double alpha_1 = U_array[2][idx_side] - U_array[1][idx_side];
                        double alpha_2 = U_array[3][idx_side] - U_array[2][idx_side];
                        double alpha_3 = U_array[4][idx_side] - U_array[3][idx_side];
                        
                        double theta_1 = fabs(alpha_1 - alpha_2)/(fabs(alpha_1) + fabs(alpha_2) + EPSILON);
                        double theta_2 = fabs(alpha_2 - alpha_3)/(fabs(alpha_2) + fabs(alpha_3) + EPSILON);
                        
                        sigma[ei][idx_side] = fmax(theta_1, theta_2);
                    }
                }
            }
        }
        
        /*
         * Compute sigma in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            sigma[ei] = variables_sigma[ei]->getPointer(1);
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(1);
            }
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1;
                     j < interior_dim_1 + 1 + num_ghosts_1;
                     j++)
                {
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_side = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0 +
                            (k + num_ghosts_2)*ghostcell_dim_0*
                                (ghostcell_dim_1 + 1);
                        
                        double alpha_1 = U_array[2][idx_side] - U_array[1][idx_side];
                        double alpha_2 = U_array[3][idx_side] - U_array[2][idx_side];
                        double alpha_3 = U_array[4][idx_side] - U_array[3][idx_side];
                        
                        double theta_1 = fabs(alpha_1 - alpha_2)/(fabs(alpha_1) + fabs(alpha_2) + EPSILON);
                        double theta_2 = fabs(alpha_2 - alpha_3)/(fabs(alpha_2) + fabs(alpha_3) + EPSILON);
                        
                        sigma[ei][idx_side] = fmax(theta_1, theta_2);
                    }
                }
            }
        }
        
        /*
         * Compute sigma in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            sigma[ei] = variables_sigma[ei]->getPointer(2);
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(2);
            }
            
            for (int k = -num_ghosts_2;
                 k < interior_dim_2 + 1 + num_ghosts_2;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_side = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0 +
                            (k + num_ghosts_2)*ghostcell_dim_0*
                                ghostcell_dim_1;
                        
                        double alpha_1 = U_array[2][idx_side] - U_array[1][idx_side];
                        double alpha_2 = U_array[3][idx_side] - U_array[2][idx_side];
                        double alpha_3 = U_array[4][idx_side] - U_array[3][idx_side];
                        
                        double theta_1 = fabs(alpha_1 - alpha_2)/(fabs(alpha_1) + fabs(alpha_2) + EPSILON);
                        double theta_2 = fabs(alpha_2 - alpha_3)/(fabs(alpha_2) + fabs(alpha_3) + EPSILON);
                        
                        sigma[ei][idx_side] = fmax(theta_1, theta_2);
                    }
                }
            }
        }
    }
}


/*
 * Compute sigma's.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::computeSigma(
    double* sigma,
    const std::vector<double*>& U_array,
    const int& idx_side)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.size()) == 6);
#endif
    
    /*
     * Compute the sigma.
     */
    
    const double alpha_1 = U_array[2][idx_side] - U_array[1][idx_side];
    const double alpha_2 = U_array[3][idx_side] - U_array[2][idx_side];
    const double alpha_3 = U_array[4][idx_side] - U_array[3][idx_side];
    
    const double theta_1 = fabs(alpha_1 - alpha_2)/(fabs(alpha_1) + fabs(alpha_2) + EPSILON);
    const double theta_2 = fabs(alpha_2 - alpha_3)/(fabs(alpha_2) + fabs(alpha_3) + EPSILON);
    
    *sigma = fmax(theta_1, theta_2);
}


/*
 * Compute beta's.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::computeBeta(
    double* beta,
    const std::vector<double*>& U_array,
    const int& idx_side)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.size()) == 6);
#endif
    
    beta[0] = 1.0/3.0*(U_array[0][idx_side]*(4.0*U_array[0][idx_side] - 19.0*U_array[1][idx_side] +
         11.0*U_array[2][idx_side]) + U_array[1][idx_side]*(25.0*U_array[1][idx_side] -
         31.0*U_array[2][idx_side]) + 10.0*U_array[2][idx_side]*U_array[2][idx_side]);
    
    beta[1] = 1.0/3.0*(U_array[1][idx_side]*(4.0*U_array[1][idx_side] - 13.0*U_array[2][idx_side] +
         5.0*U_array[3][idx_side]) + 13.0*U_array[2][idx_side]*(U_array[2][idx_side] -
         U_array[3][idx_side]) + 4.0*U_array[3][idx_side]*U_array[3][idx_side]);
    
    beta[2] = 1.0/3.0*(U_array[2][idx_side]*(10.0*U_array[2][idx_side] - 31.0*U_array[3][idx_side] +
         11.0*U_array[4][idx_side]) + U_array[3][idx_side]*(25.0*U_array[3][idx_side] -
         19.0*U_array[4][idx_side]) + 4.0*U_array[4][idx_side]*U_array[4][idx_side]);
    
    beta[3] = 1.0/232243200.0*(U_array[0][idx_side]*(525910327.0*U_array[0][idx_side] -
         4562164630.0*U_array[1][idx_side] + 7799501420.0*U_array[2][idx_side] -
         6610694540.0*U_array[3][idx_side] + 2794296070.0*U_array[4][idx_side] -
         472758974.0*U_array[5][idx_side]) + 5.0*U_array[1][idx_side]*
        (2146987907.0*U_array[1][idx_side] - 7722406988.0*U_array[2][idx_side] +
         6763559276.0*U_array[3][idx_side] - 2926461814.0*U_array[4][idx_side] +
         503766638.0*U_array[5][idx_side]) + 20.0*U_array[2][idx_side]*
        (1833221603.0*U_array[2][idx_side] - 3358664662.0*U_array[3][idx_side] +
         1495974539.0*U_array[4][idx_side] - 263126407.0*U_array[5][idx_side]) +
        20.0*U_array[3][idx_side]*(1607794163.0*U_array[3][idx_side] -
         1486026707.0*U_array[4][idx_side] + 268747951.0*U_array[5][idx_side]) +
        5.0*U_array[4][idx_side]*(1432381427.0*U_array[4][idx_side] -
         536951582.0*U_array[5][idx_side]) +
        263126407.0*U_array[5][idx_side]*U_array[5][idx_side]);
}


/*
 * Compute beta_tilde's.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::computeBetaTilde(
    double* beta_tilde,
    const std::vector<double*>& U_array,
    const int& idx_side)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.size()) == 6);
#endif
    
    beta_tilde[0] = 1.0/3.0*(U_array[5][idx_side]*(4.0*U_array[5][idx_side] - 19.0*U_array[4][idx_side] +
         11.0*U_array[3][idx_side]) + U_array[4][idx_side]*(25.0*U_array[4][idx_side] -
         31.0*U_array[3][idx_side]) + 10.0*U_array[3][idx_side]*U_array[3][idx_side]);
    
    beta_tilde[1] = 1.0/3.0*(U_array[4][idx_side]*(4.0*U_array[4][idx_side] - 13.0*U_array[3][idx_side] +
         5.0*U_array[2][idx_side]) + 13.0*U_array[3][idx_side]*(U_array[3][idx_side] -
         U_array[2][idx_side]) + 4.0*U_array[2][idx_side]*U_array[2][idx_side]);
    
    beta_tilde[2] = 1.0/3.0*(U_array[3][idx_side]*(10.0*U_array[3][idx_side] - 31.0*U_array[2][idx_side] +
         11.0*U_array[1][idx_side]) + U_array[2][idx_side]*(25.0*U_array[2][idx_side] -
         19.0*U_array[1][idx_side]) + 4.0*U_array[1][idx_side]*U_array[1][idx_side]);
    
    beta_tilde[3] = 1.0/232243200.0*(U_array[5][idx_side]*(525910327.0*U_array[5][idx_side] -
         4562164630.0*U_array[4][idx_side] + 7799501420.0*U_array[3][idx_side] -
         6610694540.0*U_array[2][idx_side] + 2794296070.0*U_array[1][idx_side] -
         472758974.0*U_array[0][idx_side]) + 5.0*U_array[4][idx_side]*
        (2146987907.0*U_array[4][idx_side] - 7722406988.0*U_array[3][idx_side] +
         6763559276.0*U_array[2][idx_side] - 2926461814.0*U_array[1][idx_side] +
         503766638.0*U_array[0][idx_side]) + 20.0*U_array[3][idx_side]*
        (1833221603.0*U_array[3][idx_side] - 3358664662.0*U_array[2][idx_side] +
         1495974539.0*U_array[1][idx_side] - 263126407.0*U_array[0][idx_side]) +
        20.0*U_array[2][idx_side]*(1607794163.0*U_array[2][idx_side] -
         1486026707.0*U_array[1][idx_side] + 268747951.0*U_array[0][idx_side]) +
        5.0*U_array[1][idx_side]*(1432381427.0*U_array[1][idx_side] -
         536951582.0*U_array[0][idx_side]) +
        263126407.0*U_array[0][idx_side]*U_array[0][idx_side]);
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::performWENOInterpolation_new(
    double* U_minus,
    double* U_plus,
    const std::vector<double*>& U_array,
    const int& idx_side)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.size()) == 6);
#endif
    
    /*
     * Perform the WENO interpolation.
     */
    
    const double& C          = d_constant_C;
    const int& p             = d_constant_p;
    const int& q             = d_constant_q;
    const double& alpha_beta = d_constant_alpha_beta;
    
    // Compute sigma.
    double sigma;
    
    #pragma forceinline
    computeSigma(&sigma, U_array, idx_side);
    
    // Compute beta's and beta_tilde's.
    double beta[4];
    double beta_tilde[4];
    
    #pragma forceinline
    computeBeta(beta, U_array, idx_side);
    
    #pragma forceinline
    computeBetaTilde(beta_tilde, U_array, idx_side);
    
    U_minus[idx_side] = beta[0] + beta[1] + beta[2] + beta[3];
    U_plus[idx_side] = beta_tilde[0] + beta_tilde[1] + beta_tilde[2] + beta_tilde[3];
}


void
ConvectiveFluxReconstructorWCNS6_Test::computeBeta(
    std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > >& variables_beta,
    const std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > >& variables_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(variables_beta.size()) == 4);
    TBOX_ASSERT(static_cast<int>(variables_beta[0].size()) == d_num_eqn);
#endif
    
    /*
     * Get the dimensions of interior box, box including ghost cells and number of ghost cells.
     */
    
    const hier::IntVector interior_dims = variables_beta[0][0]->getBox().numberCells();
    const hier::IntVector ghostcell_dims = variables_beta[0][0]->getGhostBox().numberCells();
    const hier::IntVector num_ghosts = variables_beta[0][0]->getGhostCellWidth();
    
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    for (int m = 0; m < 4; m++)
    {
        TBOX_ASSERT(static_cast<int>(variables_beta[m].size()) == d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ERROR(variables_beta[m][ei]->getBox().numberCells() == interior_dims);
            TBOX_ERROR(variables_beta[m][ei]->getGhostBox().numberCells() == ghostcell_dims);
        }
    }
    
    TBOX_ASSERT(static_cast<int>(variables_array.size()) == 6);
    
    for (int m = 0; m < 6; m++)
    {
        TBOX_ASSERT(static_cast<int>(variables_array[m].size()) == d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ERROR(variables_array[m][ei]->getBox().numberCells() == interior_dims);
            TBOX_ERROR(variables_array[m][ei]->getGhostBox().numberCells() == ghostcell_dims);
        }
    }
#endif
    
    /*
     * Declare the pointers to the beta's.
     */
    
    std::vector<double*> beta;
    beta.resize(4);
    
    /*
     * Declare the pointers to the array of variables.
     */
    
    std::vector<double*> U_array;
    U_array.resize(6);
    
    /*
     * Compute the beta's.
     */
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        const int num_ghosts_0 = num_ghosts[0];
        
        /*
         * Compute beta in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int m = 0; m < 4; m++)
            {
                beta[m] = variables_beta[m][ei]->getPointer(0);
            }
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(0);
            }
            
            #pragma ivdep
            for (int i = -num_ghosts_0;
                 i < interior_dim_0 + 1 + num_ghosts_0;
                 i++)
            {
                // Compute the linear index.
                const int idx_side = i + num_ghosts_0;
                
                beta[0][idx_side] = 1.0/3.0*(U_array[0][idx_side]*(4.0*U_array[0][idx_side] -
                    19.0*U_array[1][idx_side] + 11.0*U_array[2][idx_side]) +
                    U_array[1][idx_side]*(25.0*U_array[1][idx_side] - 31.0*U_array[2][idx_side]) +
                    10.0*U_array[2][idx_side]*U_array[2][idx_side]);
                
                beta[1][idx_side] = 1.0/3.0*(U_array[1][idx_side]*(4.0*U_array[1][idx_side] -
                    13.0*U_array[2][idx_side] + 5.0*U_array[3][idx_side]) +
                    13.0*U_array[2][idx_side]*(U_array[2][idx_side] - U_array[3][idx_side]) +
                    4.0*U_array[3][idx_side]*U_array[3][idx_side]);
                
                beta[2][idx_side] = 1.0/3.0*(U_array[2][idx_side]*(10.0*U_array[2][idx_side] -
                    31.0*U_array[3][idx_side] + 11.0*U_array[4][idx_side]) +
                    U_array[3][idx_side]*(25.0*U_array[3][idx_side] - 19.0*U_array[4][idx_side]) +
                    4.0*U_array[4][idx_side]*U_array[4][idx_side]);
                
                beta[3][idx_side] = 1.0/232243200.0*(U_array[0][idx_side]*
                    (525910327.0*U_array[0][idx_side] - 4562164630.0*U_array[1][idx_side] +
                    7799501420.0*U_array[2][idx_side] - 6610694540.0*U_array[3][idx_side] +
                    2794296070.0*U_array[4][idx_side] - 472758974.0*U_array[5][idx_side]) +
                    5.0*U_array[1][idx_side]* (2146987907.0*U_array[1][idx_side] -
                    7722406988.0*U_array[2][idx_side] + 6763559276.0*U_array[3][idx_side] -
                    2926461814.0*U_array[4][idx_side] + 503766638.0*U_array[5][idx_side]) +
                    20.0*U_array[2][idx_side]*(1833221603.0*U_array[2][idx_side] -
                    3358664662.0*U_array[3][idx_side] + 1495974539.0*U_array[4][idx_side] -
                    263126407.0*U_array[5][idx_side]) + 20.0*U_array[3][idx_side]*
                    (1607794163.0*U_array[3][idx_side] - 1486026707.0*U_array[4][idx_side] +
                    268747951.0*U_array[5][idx_side]) + 5.0*U_array[4][idx_side]*
                    (1432381427.0*U_array[4][idx_side] - 536951582.0*U_array[5][idx_side]) +
                    263126407.0*U_array[5][idx_side]*U_array[5][idx_side]);
            }
            
            #pragma ivdep
            for (int i = -num_ghosts_0;
                 i < interior_dim_0 + 1 + num_ghosts_0;
                 i++)
            {
                // Compute the linear index.
                const int idx_side = i + num_ghosts_0;
                
                double beta_test[4];
                
                #pragma forceinline
                computeBeta(beta_test, U_array, idx_side);
                
                beta[0][idx_side] = beta_test[0];
                beta[1][idx_side] = beta_test[1];
                beta[2][idx_side] = beta_test[2];
                beta[3][idx_side] = beta_test[3];
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int ghostcell_dim_0 = ghostcell_dims[0];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        
        /*
         * Compute beta in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int m = 0; m < 4; m++)
            {
                beta[m] = variables_beta[m][ei]->getPointer(0);
            }
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(0);
            }
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                #pragma ivdep
                for (int i = -num_ghosts_0;
                     i < interior_dim_1 + 1 + num_ghosts_0;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_side = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*(ghostcell_dim_0 + 1);
                    
                    beta[0][idx_side] = 1.0/3.0*(U_array[0][idx_side]*(4.0*U_array[0][idx_side] -
                        19.0*U_array[1][idx_side] + 11.0*U_array[2][idx_side]) +
                        U_array[1][idx_side]*(25.0*U_array[1][idx_side] - 31.0*U_array[2][idx_side]) +
                        10.0*U_array[2][idx_side]*U_array[2][idx_side]);
                    
                    beta[1][idx_side] = 1.0/3.0*(U_array[1][idx_side]*(4.0*U_array[1][idx_side] -
                        13.0*U_array[2][idx_side] + 5.0*U_array[3][idx_side]) +
                        13.0*U_array[2][idx_side]*(U_array[2][idx_side] - U_array[3][idx_side]) +
                        4.0*U_array[3][idx_side]*U_array[3][idx_side]);
                    
                    beta[2][idx_side] = 1.0/3.0*(U_array[2][idx_side]*(10.0*U_array[2][idx_side] -
                        31.0*U_array[3][idx_side] + 11.0*U_array[4][idx_side]) +
                        U_array[3][idx_side]*(25.0*U_array[3][idx_side] - 19.0*U_array[4][idx_side]) +
                        4.0*U_array[4][idx_side]*U_array[4][idx_side]);
                    
                    beta[3][idx_side] = 1.0/232243200.0*(U_array[0][idx_side]*
                        (525910327.0*U_array[0][idx_side] - 4562164630.0*U_array[1][idx_side] +
                        7799501420.0*U_array[2][idx_side] - 6610694540.0*U_array[3][idx_side] +
                        2794296070.0*U_array[4][idx_side] - 472758974.0*U_array[5][idx_side]) +
                        5.0*U_array[1][idx_side]* (2146987907.0*U_array[1][idx_side] -
                        7722406988.0*U_array[2][idx_side] + 6763559276.0*U_array[3][idx_side] -
                        2926461814.0*U_array[4][idx_side] + 503766638.0*U_array[5][idx_side]) +
                        20.0*U_array[2][idx_side]*(1833221603.0*U_array[2][idx_side] -
                        3358664662.0*U_array[3][idx_side] + 1495974539.0*U_array[4][idx_side] -
                        263126407.0*U_array[5][idx_side]) + 20.0*U_array[3][idx_side]*
                        (1607794163.0*U_array[3][idx_side] - 1486026707.0*U_array[4][idx_side] +
                        268747951.0*U_array[5][idx_side]) + 5.0*U_array[4][idx_side]*
                        (1432381427.0*U_array[4][idx_side] - 536951582.0*U_array[5][idx_side]) +
                        263126407.0*U_array[5][idx_side]*U_array[5][idx_side]);
                }
            }
        }
        
        /*
         * Compute beta in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int m = 0; m < 4; m++)
            {
                beta[m] = variables_beta[m][ei]->getPointer(1);
            }
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(1);
            }
            
            for (int j = -num_ghosts_1;
                 j < interior_dim_1 + 1 + num_ghosts_1;
                 j++)
            {
                #pragma ivdep
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_side = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*ghostcell_dim_0;
                    
                    beta[0][idx_side] = 1.0/3.0*(U_array[0][idx_side]*(4.0*U_array[0][idx_side] -
                        19.0*U_array[1][idx_side] + 11.0*U_array[2][idx_side]) +
                        U_array[1][idx_side]*(25.0*U_array[1][idx_side] - 31.0*U_array[2][idx_side]) +
                        10.0*U_array[2][idx_side]*U_array[2][idx_side]);
                    
                    beta[1][idx_side] = 1.0/3.0*(U_array[1][idx_side]*(4.0*U_array[1][idx_side] -
                        13.0*U_array[2][idx_side] + 5.0*U_array[3][idx_side]) +
                        13.0*U_array[2][idx_side]*(U_array[2][idx_side] - U_array[3][idx_side]) +
                        4.0*U_array[3][idx_side]*U_array[3][idx_side]);
                    
                    beta[2][idx_side] = 1.0/3.0*(U_array[2][idx_side]*(10.0*U_array[2][idx_side] -
                        31.0*U_array[3][idx_side] + 11.0*U_array[4][idx_side]) +
                        U_array[3][idx_side]*(25.0*U_array[3][idx_side] - 19.0*U_array[4][idx_side]) +
                        4.0*U_array[4][idx_side]*U_array[4][idx_side]);
                    
                    beta[3][idx_side] = 1.0/232243200.0*(U_array[0][idx_side]*
                        (525910327.0*U_array[0][idx_side] - 4562164630.0*U_array[1][idx_side] +
                        7799501420.0*U_array[2][idx_side] - 6610694540.0*U_array[3][idx_side] +
                        2794296070.0*U_array[4][idx_side] - 472758974.0*U_array[5][idx_side]) +
                        5.0*U_array[1][idx_side]* (2146987907.0*U_array[1][idx_side] -
                        7722406988.0*U_array[2][idx_side] + 6763559276.0*U_array[3][idx_side] -
                        2926461814.0*U_array[4][idx_side] + 503766638.0*U_array[5][idx_side]) +
                        20.0*U_array[2][idx_side]*(1833221603.0*U_array[2][idx_side] -
                        3358664662.0*U_array[3][idx_side] + 1495974539.0*U_array[4][idx_side] -
                        263126407.0*U_array[5][idx_side]) + 20.0*U_array[3][idx_side]*
                        (1607794163.0*U_array[3][idx_side] - 1486026707.0*U_array[4][idx_side] +
                        268747951.0*U_array[5][idx_side]) + 5.0*U_array[4][idx_side]*
                        (1432381427.0*U_array[4][idx_side] - 536951582.0*U_array[5][idx_side]) +
                        263126407.0*U_array[5][idx_side]*U_array[5][idx_side]);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int ghostcell_dim_0 = ghostcell_dims[0];
        const int ghostcell_dim_1 = ghostcell_dims[1];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        const int num_ghosts_2 = num_ghosts[2];
        
        /*
         * Compute beta in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int m = 0; m < 4; m++)
            {
                beta[m] = variables_beta[m][ei]->getPointer(0);
            }
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(0);
            }
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #pragma ivdep
                    for (int i = -num_ghosts_0;
                         i < interior_dim_1 + 1 + num_ghosts_0;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx_side = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*(ghostcell_dim_0 + 1) +
                            (k + num_ghosts_2)*(ghostcell_dim_0 + 1)*
                                ghostcell_dim_1;
                        
                        beta[0][idx_side] = 1.0/3.0*(U_array[0][idx_side]*(4.0*U_array[0][idx_side] -
                            19.0*U_array[1][idx_side] + 11.0*U_array[2][idx_side]) +
                            U_array[1][idx_side]*(25.0*U_array[1][idx_side] - 31.0*U_array[2][idx_side]) +
                            10.0*U_array[2][idx_side]*U_array[2][idx_side]);
                        
                        beta[1][idx_side] = 1.0/3.0*(U_array[1][idx_side]*(4.0*U_array[1][idx_side] -
                            13.0*U_array[2][idx_side] + 5.0*U_array[3][idx_side]) +
                            13.0*U_array[2][idx_side]*(U_array[2][idx_side] - U_array[3][idx_side]) +
                            4.0*U_array[3][idx_side]*U_array[3][idx_side]);
                        
                        beta[2][idx_side] = 1.0/3.0*(U_array[2][idx_side]*(10.0*U_array[2][idx_side] -
                            31.0*U_array[3][idx_side] + 11.0*U_array[4][idx_side]) +
                            U_array[3][idx_side]*(25.0*U_array[3][idx_side] - 19.0*U_array[4][idx_side]) +
                            4.0*U_array[4][idx_side]*U_array[4][idx_side]);
                        
                        beta[3][idx_side] = 1.0/232243200.0*(U_array[0][idx_side]*
                            (525910327.0*U_array[0][idx_side] - 4562164630.0*U_array[1][idx_side] +
                            7799501420.0*U_array[2][idx_side] - 6610694540.0*U_array[3][idx_side] +
                            2794296070.0*U_array[4][idx_side] - 472758974.0*U_array[5][idx_side]) +
                            5.0*U_array[1][idx_side]* (2146987907.0*U_array[1][idx_side] -
                            7722406988.0*U_array[2][idx_side] + 6763559276.0*U_array[3][idx_side] -
                            2926461814.0*U_array[4][idx_side] + 503766638.0*U_array[5][idx_side]) +
                            20.0*U_array[2][idx_side]*(1833221603.0*U_array[2][idx_side] -
                            3358664662.0*U_array[3][idx_side] + 1495974539.0*U_array[4][idx_side] -
                            263126407.0*U_array[5][idx_side]) + 20.0*U_array[3][idx_side]*
                            (1607794163.0*U_array[3][idx_side] - 1486026707.0*U_array[4][idx_side] +
                            268747951.0*U_array[5][idx_side]) + 5.0*U_array[4][idx_side]*
                            (1432381427.0*U_array[4][idx_side] - 536951582.0*U_array[5][idx_side]) +
                            263126407.0*U_array[5][idx_side]*U_array[5][idx_side]);
                    }
                }
            }
        }
        
        /*
         * Compute beta in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int m = 0; m < 4; m++)
            {
                beta[m] = variables_beta[m][ei]->getPointer(1);
            }
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(1);
            }
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1;
                     j < interior_dim_1 + 1 + num_ghosts_1;
                     j++)
                {
                    #pragma ivdep
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_side = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0 +
                            (k + num_ghosts_2)*ghostcell_dim_0*
                                (ghostcell_dim_1 + 1);
                        
                        beta[0][idx_side] = 1.0/3.0*(U_array[0][idx_side]*(4.0*U_array[0][idx_side] -
                            19.0*U_array[1][idx_side] + 11.0*U_array[2][idx_side]) +
                            U_array[1][idx_side]*(25.0*U_array[1][idx_side] - 31.0*U_array[2][idx_side]) +
                            10.0*U_array[2][idx_side]*U_array[2][idx_side]);
                        
                        beta[1][idx_side] = 1.0/3.0*(U_array[1][idx_side]*(4.0*U_array[1][idx_side] -
                            13.0*U_array[2][idx_side] + 5.0*U_array[3][idx_side]) +
                            13.0*U_array[2][idx_side]*(U_array[2][idx_side] - U_array[3][idx_side]) +
                            4.0*U_array[3][idx_side]*U_array[3][idx_side]);
                        
                        beta[2][idx_side] = 1.0/3.0*(U_array[2][idx_side]*(10.0*U_array[2][idx_side] -
                            31.0*U_array[3][idx_side] + 11.0*U_array[4][idx_side]) +
                            U_array[3][idx_side]*(25.0*U_array[3][idx_side] - 19.0*U_array[4][idx_side]) +
                            4.0*U_array[4][idx_side]*U_array[4][idx_side]);
                        
                        beta[3][idx_side] = 1.0/232243200.0*(U_array[0][idx_side]*
                            (525910327.0*U_array[0][idx_side] - 4562164630.0*U_array[1][idx_side] +
                            7799501420.0*U_array[2][idx_side] - 6610694540.0*U_array[3][idx_side] +
                            2794296070.0*U_array[4][idx_side] - 472758974.0*U_array[5][idx_side]) +
                            5.0*U_array[1][idx_side]* (2146987907.0*U_array[1][idx_side] -
                            7722406988.0*U_array[2][idx_side] + 6763559276.0*U_array[3][idx_side] -
                            2926461814.0*U_array[4][idx_side] + 503766638.0*U_array[5][idx_side]) +
                            20.0*U_array[2][idx_side]*(1833221603.0*U_array[2][idx_side] -
                            3358664662.0*U_array[3][idx_side] + 1495974539.0*U_array[4][idx_side] -
                            263126407.0*U_array[5][idx_side]) + 20.0*U_array[3][idx_side]*
                            (1607794163.0*U_array[3][idx_side] - 1486026707.0*U_array[4][idx_side] +
                            268747951.0*U_array[5][idx_side]) + 5.0*U_array[4][idx_side]*
                            (1432381427.0*U_array[4][idx_side] - 536951582.0*U_array[5][idx_side]) +
                            263126407.0*U_array[5][idx_side]*U_array[5][idx_side]);
                    }
                }
            }
        }
        
        /*
         * Compute beta in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int m = 0; m < 4; m++)
            {
                beta[m] = variables_beta[m][ei]->getPointer(2);
            }
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(2);
            }
            
            for (int k = -num_ghosts_2;
                 k < interior_dim_2 + 1 + num_ghosts_2;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #pragma ivdep
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_side = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0 +
                            (k + num_ghosts_2)*ghostcell_dim_0*
                                ghostcell_dim_1;
                        
                        beta[0][idx_side] = 1.0/3.0*(U_array[0][idx_side]*(4.0*U_array[0][idx_side] -
                            19.0*U_array[1][idx_side] + 11.0*U_array[2][idx_side]) +
                            U_array[1][idx_side]*(25.0*U_array[1][idx_side] - 31.0*U_array[2][idx_side]) +
                            10.0*U_array[2][idx_side]*U_array[2][idx_side]);
                        
                        beta[1][idx_side] = 1.0/3.0*(U_array[1][idx_side]*(4.0*U_array[1][idx_side] -
                            13.0*U_array[2][idx_side] + 5.0*U_array[3][idx_side]) +
                            13.0*U_array[2][idx_side]*(U_array[2][idx_side] - U_array[3][idx_side]) +
                            4.0*U_array[3][idx_side]*U_array[3][idx_side]);
                        
                        beta[2][idx_side] = 1.0/3.0*(U_array[2][idx_side]*(10.0*U_array[2][idx_side] -
                            31.0*U_array[3][idx_side] + 11.0*U_array[4][idx_side]) +
                            U_array[3][idx_side]*(25.0*U_array[3][idx_side] - 19.0*U_array[4][idx_side]) +
                            4.0*U_array[4][idx_side]*U_array[4][idx_side]);
                        
                        beta[3][idx_side] = 1.0/232243200.0*(U_array[0][idx_side]*
                            (525910327.0*U_array[0][idx_side] - 4562164630.0*U_array[1][idx_side] +
                            7799501420.0*U_array[2][idx_side] - 6610694540.0*U_array[3][idx_side] +
                            2794296070.0*U_array[4][idx_side] - 472758974.0*U_array[5][idx_side]) +
                            5.0*U_array[1][idx_side]* (2146987907.0*U_array[1][idx_side] -
                            7722406988.0*U_array[2][idx_side] + 6763559276.0*U_array[3][idx_side] -
                            2926461814.0*U_array[4][idx_side] + 503766638.0*U_array[5][idx_side]) +
                            20.0*U_array[2][idx_side]*(1833221603.0*U_array[2][idx_side] -
                            3358664662.0*U_array[3][idx_side] + 1495974539.0*U_array[4][idx_side] -
                            263126407.0*U_array[5][idx_side]) + 20.0*U_array[3][idx_side]*
                            (1607794163.0*U_array[3][idx_side] - 1486026707.0*U_array[4][idx_side] +
                            268747951.0*U_array[5][idx_side]) + 5.0*U_array[4][idx_side]*
                            (1432381427.0*U_array[4][idx_side] - 536951582.0*U_array[5][idx_side]) +
                            263126407.0*U_array[5][idx_side]*U_array[5][idx_side]);
                    }
                }
            }
        }
    }
}


void
ConvectiveFluxReconstructorWCNS6_Test::computeBetaTilde(
    std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > >& variables_beta_tilde,
    const std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > >& variables_array)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(variables_beta_tilde.size()) == 4);
    TBOX_ASSERT(static_cast<int>(variables_beta_tilde[0].size()) == d_num_eqn);
#endif
    
    /*
     * Get the dimensions of interior box and box including ghost cells.
     */
    
    const hier::IntVector interior_dims = variables_beta_tilde[0][0]->getBox().numberCells();
    const hier::IntVector ghostcell_dims = variables_beta_tilde[0][0]->getGhostBox().numberCells();
    const hier::IntVector num_ghosts = variables_beta_tilde[0][0]->getGhostCellWidth();
    
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    for (int m = 0; m < 4; m++)
    {
        TBOX_ASSERT(static_cast<int>(variables_beta_tilde[m].size()) == d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ERROR(variables_beta_tilde[m][ei]->getBox().numberCells() == interior_dims);
            TBOX_ERROR(variables_beta_tilde[m][ei]->getGhostBox().numberCells() == ghostcell_dims);
        }
    }
    
    TBOX_ASSERT(static_cast<int>(variables_array.size()) == 6);
    
    for (int m = 0; m < 6; m++)
    {
        TBOX_ASSERT(static_cast<int>(variables_array[m].size()) == d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ERROR(variables_array[m][ei]->getBox().numberCells() == interior_dims);
            TBOX_ERROR(variables_array[m][ei]->getGhostBox().numberCells() == ghostcell_dims);
        }
    }
#endif
    
    /*
     * Declare the pointers to the beta_tilde's.
     */
    
    std::vector<double*> beta_tilde;
    beta_tilde.resize(4);
    
    /*
     * Declare the pointers to the array of variables.
     */
    
    std::vector<double*> U_array;
    U_array.resize(6);
    
    /*
     * Compute the beta_tilde's.
     */
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        const int num_ghosts_0 = num_ghosts[0];
        
        /*
         * Compute beta_tilde in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int m = 0; m < 4; m++)
            {
                beta_tilde[m] = variables_beta_tilde[m][ei]->getPointer(0);
            }
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(0);
            }
            
            #pragma ivdep
            for (int i = -num_ghosts_0;
                 i < interior_dim_0 + 1 + num_ghosts_0;
                 i++)
            {
                // Compute the linear index.
                const int idx_side = i + num_ghosts_0;
                
                beta_tilde[0][idx_side] = 1.0/3.0*(U_array[5][idx_side]*(4.0*U_array[5][idx_side] -
                    19.0*U_array[4][idx_side] + 11.0*U_array[3][idx_side]) +
                    U_array[4][idx_side]*(25.0*U_array[4][idx_side] - 31.0*U_array[3][idx_side]) +
                    10.0*U_array[3][idx_side]*U_array[3][idx_side]);
                
                beta_tilde[1][idx_side] = 1.0/3.0*(U_array[4][idx_side]*(4.0*U_array[4][idx_side] -
                    13.0*U_array[3][idx_side] + 5.0*U_array[2][idx_side]) +
                    13.0*U_array[3][idx_side]*(U_array[3][idx_side] - U_array[2][idx_side]) +
                    4.0*U_array[2][idx_side]*U_array[2][idx_side]);
                
                beta_tilde[2][idx_side] = 1.0/3.0*(U_array[3][idx_side]*(10.0*U_array[3][idx_side] -
                    31.0*U_array[2][idx_side] + 11.0*U_array[1][idx_side]) +
                    U_array[2][idx_side]*(25.0*U_array[2][idx_side] - 19.0*U_array[1][idx_side]) +
                    4.0*U_array[1][idx_side]*U_array[1][idx_side]);
                
                beta_tilde[3][idx_side] = 1.0/232243200.0*(U_array[5][idx_side]*
                    (525910327.0*U_array[5][idx_side] - 4562164630.0*U_array[4][idx_side] +
                    7799501420.0*U_array[3][idx_side] - 6610694540.0*U_array[2][idx_side] +
                    2794296070.0*U_array[1][idx_side] - 472758974.0*U_array[0][idx_side]) +
                    5.0*U_array[4][idx_side]*(2146987907.0*U_array[4][idx_side] -
                    7722406988.0*U_array[3][idx_side] + 6763559276.0*U_array[2][idx_side] -
                    2926461814.0*U_array[1][idx_side] + 503766638.0*U_array[0][idx_side]) +
                    20.0*U_array[3][idx_side]*(1833221603.0*U_array[3][idx_side] -
                    3358664662.0*U_array[2][idx_side] + 1495974539.0*U_array[1][idx_side] -
                    263126407.0*U_array[0][idx_side]) + 20.0*U_array[2][idx_side]*
                    (1607794163.0*U_array[2][idx_side] - 1486026707.0*U_array[1][idx_side] +
                    268747951.0*U_array[0][idx_side]) + 5.0*U_array[1][idx_side]*
                    (1432381427.0*U_array[1][idx_side] - 536951582.0*U_array[0][idx_side]) +
                    263126407.0*U_array[0][idx_side]*U_array[0][idx_side]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int ghostcell_dim_0 = ghostcell_dims[0];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        
        /*
         * Compute beta_tilde in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int m = 0; m < 4; m++)
            {
                beta_tilde[m] = variables_beta_tilde[m][ei]->getPointer(0);
            }
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(0);
            }
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                #pragma ivdep
                for (int i = -num_ghosts_0;
                     i < interior_dim_1 + 1 + num_ghosts_0;
                     i++)
                {
                    // Compute the linear index.
                    const int idx_side = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*(ghostcell_dim_0 + 1);
                    
                    beta_tilde[0][idx_side] = 1.0/3.0*(U_array[5][idx_side]*(4.0*U_array[5][idx_side] -
                        19.0*U_array[4][idx_side] + 11.0*U_array[3][idx_side]) +
                        U_array[4][idx_side]*(25.0*U_array[4][idx_side] - 31.0*U_array[3][idx_side]) +
                        10.0*U_array[3][idx_side]*U_array[3][idx_side]);
                    
                    beta_tilde[1][idx_side] = 1.0/3.0*(U_array[4][idx_side]*(4.0*U_array[4][idx_side] -
                        13.0*U_array[3][idx_side] + 5.0*U_array[2][idx_side]) +
                        13.0*U_array[3][idx_side]*(U_array[3][idx_side] - U_array[2][idx_side]) +
                        4.0*U_array[2][idx_side]*U_array[2][idx_side]);
                    
                    beta_tilde[2][idx_side] = 1.0/3.0*(U_array[3][idx_side]*(10.0*U_array[3][idx_side] -
                        31.0*U_array[2][idx_side] + 11.0*U_array[1][idx_side]) +
                        U_array[2][idx_side]*(25.0*U_array[2][idx_side] - 19.0*U_array[1][idx_side]) +
                        4.0*U_array[1][idx_side]*U_array[1][idx_side]);
                    
                    beta_tilde[3][idx_side] = 1.0/232243200.0*(U_array[5][idx_side]*
                        (525910327.0*U_array[5][idx_side] - 4562164630.0*U_array[4][idx_side] +
                        7799501420.0*U_array[3][idx_side] - 6610694540.0*U_array[2][idx_side] +
                        2794296070.0*U_array[1][idx_side] - 472758974.0*U_array[0][idx_side]) +
                        5.0*U_array[4][idx_side]*(2146987907.0*U_array[4][idx_side] -
                        7722406988.0*U_array[3][idx_side] + 6763559276.0*U_array[2][idx_side] -
                        2926461814.0*U_array[1][idx_side] + 503766638.0*U_array[0][idx_side]) +
                        20.0*U_array[3][idx_side]*(1833221603.0*U_array[3][idx_side] -
                        3358664662.0*U_array[2][idx_side] + 1495974539.0*U_array[1][idx_side] -
                        263126407.0*U_array[0][idx_side]) + 20.0*U_array[2][idx_side]*
                        (1607794163.0*U_array[2][idx_side] - 1486026707.0*U_array[1][idx_side] +
                        268747951.0*U_array[0][idx_side]) + 5.0*U_array[1][idx_side]*
                        (1432381427.0*U_array[1][idx_side] - 536951582.0*U_array[0][idx_side]) +
                        263126407.0*U_array[0][idx_side]*U_array[0][idx_side]);
                }
            }
        }
        
        /*
         * Compute beta_tilde in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int m = 0; m < 4; m++)
            {
                beta_tilde[m] = variables_beta_tilde[m][ei]->getPointer(1);
            }
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(1);
            }
            
            for (int j = -num_ghosts_1;
                 j < interior_dim_1 + 1 + num_ghosts_1;
                 j++)
            {
                #pragma ivdep
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index.
                    const int idx_side = (i + num_ghosts_0) +
                        (j + num_ghosts_1)*ghostcell_dim_0;
                    
                    beta_tilde[0][idx_side] = 1.0/3.0*(U_array[5][idx_side]*(4.0*U_array[5][idx_side] -
                        19.0*U_array[4][idx_side] + 11.0*U_array[3][idx_side]) +
                        U_array[4][idx_side]*(25.0*U_array[4][idx_side] - 31.0*U_array[3][idx_side]) +
                        10.0*U_array[3][idx_side]*U_array[3][idx_side]);
                    
                    beta_tilde[1][idx_side] = 1.0/3.0*(U_array[4][idx_side]*(4.0*U_array[4][idx_side] -
                        13.0*U_array[3][idx_side] + 5.0*U_array[2][idx_side]) +
                        13.0*U_array[3][idx_side]*(U_array[3][idx_side] - U_array[2][idx_side]) +
                        4.0*U_array[2][idx_side]*U_array[2][idx_side]);
                    
                    beta_tilde[2][idx_side] = 1.0/3.0*(U_array[3][idx_side]*(10.0*U_array[3][idx_side] -
                        31.0*U_array[2][idx_side] + 11.0*U_array[1][idx_side]) +
                        U_array[2][idx_side]*(25.0*U_array[2][idx_side] - 19.0*U_array[1][idx_side]) +
                        4.0*U_array[1][idx_side]*U_array[1][idx_side]);
                    
                    beta_tilde[3][idx_side] = 1.0/232243200.0*(U_array[5][idx_side]*
                        (525910327.0*U_array[5][idx_side] - 4562164630.0*U_array[4][idx_side] +
                        7799501420.0*U_array[3][idx_side] - 6610694540.0*U_array[2][idx_side] +
                        2794296070.0*U_array[1][idx_side] - 472758974.0*U_array[0][idx_side]) +
                        5.0*U_array[4][idx_side]*(2146987907.0*U_array[4][idx_side] -
                        7722406988.0*U_array[3][idx_side] + 6763559276.0*U_array[2][idx_side] -
                        2926461814.0*U_array[1][idx_side] + 503766638.0*U_array[0][idx_side]) +
                        20.0*U_array[3][idx_side]*(1833221603.0*U_array[3][idx_side] -
                        3358664662.0*U_array[2][idx_side] + 1495974539.0*U_array[1][idx_side] -
                        263126407.0*U_array[0][idx_side]) + 20.0*U_array[2][idx_side]*
                        (1607794163.0*U_array[2][idx_side] - 1486026707.0*U_array[1][idx_side] +
                        268747951.0*U_array[0][idx_side]) + 5.0*U_array[1][idx_side]*
                        (1432381427.0*U_array[1][idx_side] - 536951582.0*U_array[0][idx_side]) +
                        263126407.0*U_array[0][idx_side]*U_array[0][idx_side]);
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int ghostcell_dim_0 = ghostcell_dims[0];
        const int ghostcell_dim_1 = ghostcell_dims[1];
        
        const int num_ghosts_0 = num_ghosts[0];
        const int num_ghosts_1 = num_ghosts[1];
        const int num_ghosts_2 = num_ghosts[2];
        
        /*
         * Compute beta_tilde in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int m = 0; m < 4; m++)
            {
                beta_tilde[m] = variables_beta_tilde[m][ei]->getPointer(0);
            }
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(0);
            }
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #pragma ivdep
                    for (int i = -num_ghosts_0;
                         i < interior_dim_1 + 1 + num_ghosts_0;
                         i++)
                    {
                        // Compute the linear index.
                        const int idx_side = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*(ghostcell_dim_0 + 1) +
                            (k + num_ghosts_2)*(ghostcell_dim_0 + 1)*
                                ghostcell_dim_1;
                        
                        beta_tilde[0][idx_side] = 1.0/3.0*(U_array[5][idx_side]*(4.0*U_array[5][idx_side] -
                            19.0*U_array[4][idx_side] + 11.0*U_array[3][idx_side]) +
                            U_array[4][idx_side]*(25.0*U_array[4][idx_side] - 31.0*U_array[3][idx_side]) +
                            10.0*U_array[3][idx_side]*U_array[3][idx_side]);
                        
                        beta_tilde[1][idx_side] = 1.0/3.0*(U_array[4][idx_side]*(4.0*U_array[4][idx_side] -
                            13.0*U_array[3][idx_side] + 5.0*U_array[2][idx_side]) +
                            13.0*U_array[3][idx_side]*(U_array[3][idx_side] - U_array[2][idx_side]) +
                            4.0*U_array[2][idx_side]*U_array[2][idx_side]);
                        
                        beta_tilde[2][idx_side] = 1.0/3.0*(U_array[3][idx_side]*(10.0*U_array[3][idx_side] -
                            31.0*U_array[2][idx_side] + 11.0*U_array[1][idx_side]) +
                            U_array[2][idx_side]*(25.0*U_array[2][idx_side] - 19.0*U_array[1][idx_side]) +
                            4.0*U_array[1][idx_side]*U_array[1][idx_side]);
                        
                        beta_tilde[3][idx_side] = 1.0/232243200.0*(U_array[5][idx_side]*
                            (525910327.0*U_array[5][idx_side] - 4562164630.0*U_array[4][idx_side] +
                            7799501420.0*U_array[3][idx_side] - 6610694540.0*U_array[2][idx_side] +
                            2794296070.0*U_array[1][idx_side] - 472758974.0*U_array[0][idx_side]) +
                            5.0*U_array[4][idx_side]*(2146987907.0*U_array[4][idx_side] -
                            7722406988.0*U_array[3][idx_side] + 6763559276.0*U_array[2][idx_side] -
                            2926461814.0*U_array[1][idx_side] + 503766638.0*U_array[0][idx_side]) +
                            20.0*U_array[3][idx_side]*(1833221603.0*U_array[3][idx_side] -
                            3358664662.0*U_array[2][idx_side] + 1495974539.0*U_array[1][idx_side] -
                            263126407.0*U_array[0][idx_side]) + 20.0*U_array[2][idx_side]*
                            (1607794163.0*U_array[2][idx_side] - 1486026707.0*U_array[1][idx_side] +
                            268747951.0*U_array[0][idx_side]) + 5.0*U_array[1][idx_side]*
                            (1432381427.0*U_array[1][idx_side] - 536951582.0*U_array[0][idx_side]) +
                            263126407.0*U_array[0][idx_side]*U_array[0][idx_side]);
                    }
                }
            }
        }
        
        /*
         * Compute beta_tilde in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int m = 0; m < 4; m++)
            {
                beta_tilde[m] = variables_beta_tilde[m][ei]->getPointer(1);
            }
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(1);
            }
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -num_ghosts_1;
                     j < interior_dim_1 + 1 + num_ghosts_1;
                     j++)
                {
                    #pragma ivdep
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_side = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0 +
                            (k + num_ghosts_2)*ghostcell_dim_0*
                                (ghostcell_dim_1 + 1);
                        
                        beta_tilde[0][idx_side] = 1.0/3.0*(U_array[5][idx_side]*(4.0*U_array[5][idx_side] -
                            19.0*U_array[4][idx_side] + 11.0*U_array[3][idx_side]) +
                            U_array[4][idx_side]*(25.0*U_array[4][idx_side] - 31.0*U_array[3][idx_side]) +
                            10.0*U_array[3][idx_side]*U_array[3][idx_side]);
                        
                        beta_tilde[1][idx_side] = 1.0/3.0*(U_array[4][idx_side]*(4.0*U_array[4][idx_side] -
                            13.0*U_array[3][idx_side] + 5.0*U_array[2][idx_side]) +
                            13.0*U_array[3][idx_side]*(U_array[3][idx_side] - U_array[2][idx_side]) +
                            4.0*U_array[2][idx_side]*U_array[2][idx_side]);
                        
                        beta_tilde[2][idx_side] = 1.0/3.0*(U_array[3][idx_side]*(10.0*U_array[3][idx_side] -
                            31.0*U_array[2][idx_side] + 11.0*U_array[1][idx_side]) +
                            U_array[2][idx_side]*(25.0*U_array[2][idx_side] - 19.0*U_array[1][idx_side]) +
                            4.0*U_array[1][idx_side]*U_array[1][idx_side]);
                        
                        beta_tilde[3][idx_side] = 1.0/232243200.0*(U_array[5][idx_side]*
                            (525910327.0*U_array[5][idx_side] - 4562164630.0*U_array[4][idx_side] +
                            7799501420.0*U_array[3][idx_side] - 6610694540.0*U_array[2][idx_side] +
                            2794296070.0*U_array[1][idx_side] - 472758974.0*U_array[0][idx_side]) +
                            5.0*U_array[4][idx_side]*(2146987907.0*U_array[4][idx_side] -
                            7722406988.0*U_array[3][idx_side] + 6763559276.0*U_array[2][idx_side] -
                            2926461814.0*U_array[1][idx_side] + 503766638.0*U_array[0][idx_side]) +
                            20.0*U_array[3][idx_side]*(1833221603.0*U_array[3][idx_side] -
                            3358664662.0*U_array[2][idx_side] + 1495974539.0*U_array[1][idx_side] -
                            263126407.0*U_array[0][idx_side]) + 20.0*U_array[2][idx_side]*
                            (1607794163.0*U_array[2][idx_side] - 1486026707.0*U_array[1][idx_side] +
                            268747951.0*U_array[0][idx_side]) + 5.0*U_array[1][idx_side]*
                            (1432381427.0*U_array[1][idx_side] - 536951582.0*U_array[0][idx_side]) +
                            263126407.0*U_array[0][idx_side]*U_array[0][idx_side]);
                    }
                }
            }
        }
        
        /*
         * Compute beta_tilde in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int m = 0; m < 4; m++)
            {
                beta_tilde[m] = variables_beta_tilde[m][ei]->getPointer(2);
            }
            
            for (int m = 0; m < 6; m++)
            {
                U_array[m] = variables_array[m][ei]->getPointer(2);
            }
            
            for (int k = -num_ghosts_2;
                 k < interior_dim_2 + 1 + num_ghosts_2;
                 k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #pragma ivdep
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index.
                        const int idx_side = (i + num_ghosts_0) +
                            (j + num_ghosts_1)*ghostcell_dim_0 +
                            (k + num_ghosts_2)*ghostcell_dim_0*
                                ghostcell_dim_1;
                        
                        beta_tilde[0][idx_side] = 1.0/3.0*(U_array[5][idx_side]*(4.0*U_array[5][idx_side] -
                            19.0*U_array[4][idx_side] + 11.0*U_array[3][idx_side]) +
                            U_array[4][idx_side]*(25.0*U_array[4][idx_side] - 31.0*U_array[3][idx_side]) +
                            10.0*U_array[3][idx_side]*U_array[3][idx_side]);
                        
                        beta_tilde[1][idx_side] = 1.0/3.0*(U_array[4][idx_side]*(4.0*U_array[4][idx_side] -
                            13.0*U_array[3][idx_side] + 5.0*U_array[2][idx_side]) +
                            13.0*U_array[3][idx_side]*(U_array[3][idx_side] - U_array[2][idx_side]) +
                            4.0*U_array[2][idx_side]*U_array[2][idx_side]);
                        
                        beta_tilde[2][idx_side] = 1.0/3.0*(U_array[3][idx_side]*(10.0*U_array[3][idx_side] -
                            31.0*U_array[2][idx_side] + 11.0*U_array[1][idx_side]) +
                            U_array[2][idx_side]*(25.0*U_array[2][idx_side] - 19.0*U_array[1][idx_side]) +
                            4.0*U_array[1][idx_side]*U_array[1][idx_side]);
                        
                        beta_tilde[3][idx_side] = 1.0/232243200.0*(U_array[5][idx_side]*
                            (525910327.0*U_array[5][idx_side] - 4562164630.0*U_array[4][idx_side] +
                            7799501420.0*U_array[3][idx_side] - 6610694540.0*U_array[2][idx_side] +
                            2794296070.0*U_array[1][idx_side] - 472758974.0*U_array[0][idx_side]) +
                            5.0*U_array[4][idx_side]*(2146987907.0*U_array[4][idx_side] -
                            7722406988.0*U_array[3][idx_side] + 6763559276.0*U_array[2][idx_side] -
                            2926461814.0*U_array[1][idx_side] + 503766638.0*U_array[0][idx_side]) +
                            20.0*U_array[3][idx_side]*(1833221603.0*U_array[3][idx_side] -
                            3358664662.0*U_array[2][idx_side] + 1495974539.0*U_array[1][idx_side] -
                            263126407.0*U_array[0][idx_side]) + 20.0*U_array[2][idx_side]*
                            (1607794163.0*U_array[2][idx_side] - 1486026707.0*U_array[1][idx_side] +
                            268747951.0*U_array[0][idx_side]) + 5.0*U_array[1][idx_side]*
                            (1432381427.0*U_array[1][idx_side] - 536951582.0*U_array[0][idx_side]) +
                            263126407.0*U_array[0][idx_side]*U_array[0][idx_side]);
                    }
                }
            }
        }
    }
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::performWENOInterpolation(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_minus,
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_plus,
    const std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > >& variables_array)
{
    
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::performWENOInterpolation(
    std::vector<double>& U_minus,
    std::vector<double>& U_plus,
    const boost::multi_array<const double*, 2>& U_array,
    const hier::Index& cell_index_minus,
    const hier::Index& cell_index_plus,
    const DIRECTION::TYPE& direction)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.shape()[0]) == 6);
    TBOX_ASSERT(static_cast<int>(U_array.shape()[1]) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(U_minus.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(U_plus.size()) == d_num_eqn);
#endif
    
    /*
     * Compute the projection matrix.
     * Transform the physical variables into the characteristic variables.
     */
    
    d_flow_model->computeLocalFaceProjectionMatrixOfPrimitiveVariables(
        R_inv_intercell,
        cell_index_minus,
        cell_index_plus,
        direction);
    
    projectPhysicalVariablesToCharacteristicFields(W_array, U_array, R_inv_intercell);
    
    /*
     * Perform the WENO interpolation.
     */
    
    const double& C = d_constant_C;
    const int& p = d_constant_p;
    const int& q = d_constant_q;
    const double& alpha_beta = d_constant_alpha_beta;
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        boost::multi_array_ref<double, 2>::const_array_view<1>::type W_array_ei =
            W_array[boost::indices[boost::multi_array_ref<double, 2>::index_range()][ei]];
        
        // Compute sigma.
        double sigma;
        computeSigma(sigma, W_array_ei);
        
        // Compute beta's.
        computeBeta(beta, W_array_ei);
        computeBetaTilde(beta_tilde, W_array_ei);
        
        /*
         * Compute W_minus of the current characteristic variable.
         */
        
        // Compute the reference smoothness indicators tau_6.
        const double beta_avg = 1.0/8*(beta[0] + beta[2] + 6*beta[1]);
        const double tau_6 = fabs(beta[3] - beta_avg);
        
        if(fabs(tau_6/(beta_avg + EPSILON)) > alpha_beta)
        {
            /*
             * Compute the weights alpha_upwind.
             */
            
            double alpha_upwind[4];
            double alpha_upwind_sum = 0.0;
            
            // Define linear weights d.
            double d[4];
            d[0] = 1.0/16.0;
            d[1] = 5.0/8.0;
            d[2] = 5.0/16.0;
            
            const double tau_5 = fabs(beta[0] - beta[2]);
            
            for (int r = 0; r < 3; r++)
            {
                // Compute the weights alpha_upwind.
                alpha_upwind[r] = d[r]*(1.0 + pow(tau_5/(beta[r] + EPSILON), p));
                
                // Sum up the weights alpha_upwind.
                alpha_upwind_sum += alpha_upwind[r];
            }
            alpha_upwind[3] = 0.0;
            
            /*
             * Compute the weights alpha_central.
             */
            
            double alpha_central[4];
            double alpha_central_sum = 0.0;
            
            // Define linear weights d.
            d[0] = 1.0/32.0;
            d[1] = 15.0/32.0;
            d[2] = 15.0/32.0;
            d[3] = 1.0/32.0;
            
            for (int r = 0; r < 4; r++)
            {
                // Compute the weights alpha.
                alpha_central[r] = d[r]*(C + pow(tau_6/(beta[r] + EPSILON), q));
                
                // Sum up the weights alpha.
                alpha_central_sum += alpha_central[r];
            }
            
            // Compute the W_minus.
            W_minus[ei] = 0.0;
            
            for (int r = 0; r < 4; r++)
            {
                // Compute the linear interpolated value.
                double W_minus_r = 0.0;
                for (int m = r; m < 3 + r; m++)
                {
                    W_minus_r += d_weights_c[r][m - r]*W_array[m][ei];
                }
                
                // Compute omega.
                const double omega_upwind = alpha_upwind[r]/alpha_upwind_sum;
                const double omega_central = alpha_central[r]/alpha_central_sum;
                const double omega = sigma*omega_upwind + (1.0 - sigma)*omega_central;
                
                // Compute the nonlinear interpolated value.
                W_minus[ei] += omega*W_minus_r;
            }
        }
        else
        {
            // Define linear weights d.
            double d[4];
            d[0] = 1.0/32.0;
            d[1] = 15.0/32.0;
            d[2] = 15.0/32.0;
            d[3] = 1.0/32.0;
            
            /*
             * Compute the weights alpha.
             */
            
            double alpha[4];
            double alpha_sum = 0.0;
            
            for (int r = 0; r < 4; r++)
            {
                // Compute the weights alpha.
                alpha[r] = d[r]*(C + pow(tau_6/(beta[r] + EPSILON), q));
                
                // Sum up the weights alpha.
                alpha_sum += alpha[r];
            }
            
            // Compute the W_minus.
            W_minus[ei] = 0.0;
            
            for (int r = 0; r < 4; r++)
            {
                // Compute the linear interpolated value.
                double W_minus_r = 0.0;
                for (int m = r; m < 3 + r; m++)
                {
                    W_minus_r += d_weights_c[r][m - r]*W_array[m][ei];
                }
                
                // Compute omega.
                const double omega = alpha[r]/alpha_sum;
                
                // Compute the nonlinear interpolated value.
                W_minus[ei] += omega*W_minus_r;
            }
        }
        
        /*
         * Compute W_plus of the current characteristic variable.
         */
        
        // Compute the reference smoothness indicators tau_6_tilde.
        const double beta_tilde_avg =  1.0/8*(beta_tilde[0] + beta_tilde[2] + 6*beta_tilde[1]);
        const double tau_6_tilde = fabs(beta_tilde[3] - beta_tilde_avg);
        
        if (fabs(tau_6_tilde/(beta_tilde_avg + EPSILON)) > alpha_beta)
        {
            /*
             * Compute the weights alpha_upwind_tilde.
             */
            
            double alpha_upwind_tilde[4];
            double alpha_upwind_tilde_sum = 0.0;
            
            // Define linear weights d.
            double d[4];
            d[0] = 1.0/16.0;
            d[1] = 5.0/8.0;
            d[2] = 5.0/16.0;
            
            const double tau_5_tilde = fabs(beta_tilde[0] - beta_tilde[2]);
            
            for (int r = 0; r < 3; r++)
            {
                // Compute the weights alpha_upwind_tilde.
                alpha_upwind_tilde[r] = d[r]*(1.0 + pow(tau_5_tilde/(beta_tilde[r] + EPSILON), p));
                
                // Sum up the weights alpha_upwind_tilde.
                alpha_upwind_tilde_sum += alpha_upwind_tilde[r];
            }
            alpha_upwind_tilde[3] = 0.0;
            
            /*
             * Compute the weights alpha_central_tilde.
             */
            
            double alpha_central_tilde[4];
            double alpha_central_tilde_sum = 0.0;
            
            // Define linear weights d.
            d[0] = 1.0/32.0;
            d[1] = 15.0/32.0;
            d[2] = 15.0/32.0;
            d[3] = 1.0/32.0;
            
            for (int r = 0; r < 4; r++)
            {
                // Compute the weights alpha_tilde.
                alpha_central_tilde[r] = d[r]*(C + pow(tau_6_tilde/(beta_tilde[r] + EPSILON), q));
                
                // Sum up the weights alpha.
                alpha_central_tilde_sum += alpha_central_tilde[r];
            }
            
            // Compute the W_plus.
            W_plus[ei] = 0.0;
            
            for (int r = 0; r < 4; r++)
            {
                // Compute the linear interpolated value.
                double W_plus_r = 0.0;
                for (int m = r; m < 3 + r; m++)
                {
                    W_plus_r += d_weights_c[r][m - r]*W_array[6 - m - 1][ei];
                }
                
                // Compute omega_tilde;
                const double omega_upwind_tilde = alpha_upwind_tilde[r]/alpha_upwind_tilde_sum;
                const double omega_central_tilde = alpha_central_tilde[r]/alpha_central_tilde_sum;
                const double omega_tilde = sigma*omega_upwind_tilde + (1.0 - sigma)*omega_central_tilde;
                
                // Compute the nonlinear interpolated value.
                W_plus[ei] += omega_tilde*W_plus_r;
            }
        }
        else
        {
            // Define linear weights d.
            double d[4];
            d[0] = 1.0/32.0;
            d[1] = 15.0/32.0;
            d[2] = 15.0/32.0;
            d[3] = 1.0/32.0;
            
            /*
             * Compute the weights alpha_tilde.
             */
            
            double alpha_tilde[4];
            double alpha_tilde_sum = 0.0;
            
            for (int r = 0; r < 4; r++)
            {
                // Compute the weights alpha_tilde.
                alpha_tilde[r] = d[r]*(C + pow(tau_6_tilde/(beta_tilde[r] + EPSILON), q));
                
                // Sum up the weights alpha.
                alpha_tilde_sum += alpha_tilde[r];
            }
            
            // Compute the W_plus.
            W_plus[ei] = 0.0;
            
            for (int r = 0; r < 4; r++)
            {
                // Compute the linear interpolated value.
                double W_plus_r = 0.0;
                for (int m = r; m < 3 + r; m++)
                {
                    W_plus_r += d_weights_c[r][m - r]*W_array[6 - m - 1][ei];
                }
                
                // Compute omega_tilde;
                const double omega_tilde = alpha_tilde[r]/alpha_tilde_sum;
                
                // Compute the nonlinear interpolated value.
                W_plus[ei] += omega_tilde*W_plus_r;
            }
        }
    }
    
    /*
     * Compute the inverse of projection matrix.
     * Transform the characteristic variables back to physcial variables.
     */
    
    d_flow_model->computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables(
        R_intercell,
        cell_index_minus,
        cell_index_plus,
        direction);
    
    projectCharacteristicVariablesToPhysicalFields(U_minus, W_minus, R_intercell);
    projectCharacteristicVariablesToPhysicalFields(U_plus, W_plus, R_intercell);
    
}
