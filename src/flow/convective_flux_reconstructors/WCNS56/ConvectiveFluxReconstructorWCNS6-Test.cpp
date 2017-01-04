#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS6-Test.hpp"

#define EPSILON 1e-40

/*
 * Timers interspersed throughout the class.
 */

boost::shared_ptr<tbox::Timer> ConvectiveFluxReconstructorWCNS6_Test::t_characteristic_decomposition;
boost::shared_ptr<tbox::Timer> ConvectiveFluxReconstructorWCNS6_Test::t_WENO_interpolation;
boost::shared_ptr<tbox::Timer> ConvectiveFluxReconstructorWCNS6_Test::t_Riemann_solver;
boost::shared_ptr<tbox::Timer> ConvectiveFluxReconstructorWCNS6_Test::t_reconstruct_flux;
boost::shared_ptr<tbox::Timer> ConvectiveFluxReconstructorWCNS6_Test::t_compute_source;

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
    
    d_constant_alpha_tau = d_convective_flux_reconstructor_db->getDoubleWithDefault("constant_alpha_tau", 35.0);
    d_constant_alpha_tau = d_convective_flux_reconstructor_db->getDoubleWithDefault("d_constant_alpha_tau", d_constant_alpha_tau);
    
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
    
    t_compute_source = tbox::TimerManager::getManager()->
        getTimer("ConvectiveFluxReconstructorWCNS6_Test::t_compute_source");
    
}


ConvectiveFluxReconstructorWCNS6_Test::~ConvectiveFluxReconstructorWCNS6_Test()
{
    t_characteristic_decomposition.reset();
    t_WENO_interpolation.reset();
    t_Riemann_solver.reset();
    t_reconstruct_flux.reset();
    t_compute_source.reset();
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
    os << "d_constant_alpha_tau = "
       << d_constant_alpha_tau
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
    restart_db->putDouble("d_constant_alpha_tau", d_constant_alpha_tau);
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
    
    boost::shared_ptr<pdat::SideData<double> > convective_flux_midpoint(
        new pdat::SideData<double>(interior_box, d_num_eqn, hier::IntVector::getOne(d_dim)));
    
    boost::shared_ptr<pdat::CellData<double> > vorticity_magnitude(
        new pdat::CellData<double>(interior_box, 1, d_num_conv_ghosts));
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimension.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
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
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        
        const int num_subghosts_0_velocity = num_subghosts_velocity[0];
        const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
        
        double* u = velocity->getPointer(0);
        std::vector<double*> F_node_x;
        F_node_x.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
        }
        std::vector<double*> F_midpoint_x;
        F_midpoint_x.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
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
                
                count_eqn++;
            }
        }
        
        /*
         * Declare temporary data containers for WENO interpolation.
         */
        
        std::vector<boost::shared_ptr<pdat::SideData<double> > > projection_variables;
        
        std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > > characteristic_variables;
        
        std::vector<boost::shared_ptr<pdat::SideData<double> > > characteristic_variables_minus;
        std::vector<boost::shared_ptr<pdat::SideData<double> > > characteristic_variables_plus;
        
        std::vector<boost::shared_ptr<pdat::SideData<double> > > primitive_variables_minus;
        std::vector<boost::shared_ptr<pdat::SideData<double> > > primitive_variables_plus;
        
        /*
         * Initialize temporary data containers for WENO interpolation.
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
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> W_array;
            W_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                W_array.push_back(characteristic_variables[m][ei]->getPointer(0));
            }
            
            double* W_L = characteristic_variables_minus[ei]->getPointer(0);
            double* W_R = characteristic_variables_plus[ei]->getPointer(0);
            
            #ifdef __INTEL_COMPILER
            #pragma ivdep
            #endif
            for (int i = -1; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear index of the mid-point.
                const int idx_midpoint_x = i + 1;
                
                #ifdef __INTEL_COMPILER
                #pragma forceinline
                #endif
                performWENOInterpolation_minus(
                    W_L,
                    W_array,
                    idx_midpoint_x);
            }
            
            #ifdef __INTEL_COMPILER
            #pragma ivdep
            #endif
            for (int i = -1; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear index of the mid-point.
                const int idx_midpoint_x = i + 1;
                
                #ifdef __INTEL_COMPILER
                #pragma forceinline
                #endif
                performWENOInterpolation_plus(
                    W_R,
                    W_array,
                    idx_midpoint_x);
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
        
        {
            std::vector<double*> V_L;
            std::vector<double*> V_R;
            V_L.reserve(d_num_eqn);
            V_R.reserve(d_num_eqn);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_L.push_back(primitive_variables_minus[ei]->getPointer(0));
                V_R.push_back(primitive_variables_plus[ei]->getPointer(0));
            }
            
            // Declare and initialize containers to store the references to the
            // WENO interpolated values.
            std::vector<boost::reference_wrapper<double> > V_L_ref;
            std::vector<boost::reference_wrapper<double> > V_R_ref;
            V_L_ref.reserve(d_num_eqn);
            V_R_ref.reserve(d_num_eqn);
            
            // Declare container to store the references to the mid-point flux.
            std::vector<boost::reference_wrapper<double> > F_midpoint_x_ref;
            F_midpoint_x_ref.reserve(d_num_eqn);
            
            // Declare container to store the references to the mid-point velocity.
            std::vector<boost::reference_wrapper<double> > vel_midpoint_x_ref;
            vel_midpoint_x_ref.reserve(d_dim.getValue());
            
            for (int i = -1; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear index of the side.
                const int idx_midpoint_x = i + 1;
                
                // Initialize container that stores the references to the mid-point flux.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    F_midpoint_x_ref.push_back(boost::ref(F_midpoint_x[ei][idx_midpoint_x]));
                }
                
                // Initialize container that stores the references to the mid-point velocity.
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    vel_midpoint_x_ref.push_back(
                        boost::ref(velocity_intercell->getPointer(0, di)[idx_midpoint_x]));
                }
                
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    V_L_ref.push_back(boost::ref(V_L[ei][idx_midpoint_x]));
                    V_R_ref.push_back(boost::ref(V_R[ei][idx_midpoint_x]));
                }
                
                d_flow_model->
                    computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                        F_midpoint_x_ref,
                        vel_midpoint_x_ref,
                        V_L_ref,
                        V_R_ref,
                        DIRECTION::X_DIRECTION,
                        RIEMANN_SOLVER::HLLC);
                
                F_midpoint_x_ref.clear();
                vel_midpoint_x_ref.clear();
                V_L_ref.clear();
                V_R_ref.clear();
            }
        }
        
        t_Riemann_solver->stop();
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        t_reconstruct_flux->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            double* F_face_x = convective_flux->getPointer(0, ei);
            
            #ifdef __INTEL_COMPILER
            #pragma ivdep
            #endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                const int idx_midpoint_x = i + 1;
                const int idx_midpoint_x_L = i;
                const int idx_midpoint_x_R = i + 2;
                const int idx_node_L = i - 1 + num_subghosts_0_convective_flux_x;
                const int idx_node_R = i + num_subghosts_0_convective_flux_x;
                
                F_face_x[idx_face_x] = dt*(1.0/30.0*(F_midpoint_x[ei][idx_midpoint_x_R] +
                     F_midpoint_x[ei][idx_midpoint_x_L]) -
                    3.0/10.0*(F_node_x[ei][idx_node_R] +
                     F_node_x[ei][idx_node_L]) +
                    23.0/15.0*F_midpoint_x[ei][idx_midpoint_x]);
            }
        }
        
        t_reconstruct_flux->stop();
        
        /*
         * Compute the source.
         */
        
        t_compute_source->start();
        
        const std::vector<EQN_FORM::TYPE> eqn_form = d_flow_model->getEquationsForm();
        
        double* u_midpoint_x = velocity_intercell->getPointer(0, 0);
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            if (eqn_form[ei] == EQN_FORM::ADVECTIVE)
            {
                double* S = source->getPointer(ei);
                
                const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                
                #ifdef __INTEL_COMPILER
                #pragma ivdep
                #endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices. 
                    const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                    
                    const int idx_cell_wghost_x_L = i - 1 + num_subghosts_0_velocity;
                    
                    const int idx_cell_wghost_x_R = i + 1 +num_subghosts_0_velocity;
                    
                    const int idx_cell_nghost = i;
                    
                    const int idx_midpoint_x_LL = i;
                    
                    const int idx_midpoint_x_L = i + 1;
                    
                    const int idx_midpoint_x_R = i + 2;
                    
                    const int idx_midpoint_x_RR = i + 3;
                    
                    S[idx_cell_nghost] += dt*Q[ei][idx_cell_wghost]*(
                        (3.0/2.0*(u_midpoint_x[idx_midpoint_x_R] - u_midpoint_x[idx_midpoint_x_L]) -
                         3.0/10.0*(u[idx_cell_wghost_x_R] - u[idx_cell_wghost_x_L]) +
                         1.0/30.0*(u_midpoint_x[idx_midpoint_x_RR] - u_midpoint_x[idx_midpoint_x_LL]))/dx[0]);
                }
            }
        }
        
        t_compute_source->stop();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_conv_ghosts_0 = d_num_conv_ghosts[0];
        const int num_conv_ghosts_1 = d_num_conv_ghosts[1];
        
        const int conv_ghostcell_dim_0 = conv_ghostcell_dims[0];
        
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
        
        const int num_subghosts_0_velocity = num_subghosts_velocity[0];
        const int num_subghosts_1_velocity = num_subghosts_velocity[1];
        const int subghostcell_dim_0_velocity = subghostcell_dims_velocity[0];
        
        const int num_subghosts_0_dilatation = num_subghosts_dilatation[0];
        const int num_subghosts_1_dilatation = num_subghosts_dilatation[1];
        const int subghostcell_dim_0_dilatation = subghostcell_dims_dilatation[0];
        
        const int num_subghosts_0_vorticity = num_subghosts_vorticity[0];
        const int num_subghosts_1_vorticity = num_subghosts_vorticity[1];
        const int subghostcell_dim_0_vorticity = subghostcell_dims_vorticity[0];
        
        const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
        const int num_subghosts_1_convective_flux_x = num_subghosts_convective_flux_x[1];
        const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
        
        const int num_subghosts_0_convective_flux_y = num_subghosts_convective_flux_y[0];
        const int num_subghosts_1_convective_flux_y = num_subghosts_convective_flux_y[1];
        const int subghostcell_dim_0_convective_flux_y = subghostcell_dims_convective_flux_y[0];
        
        double* u     = velocity->getPointer(0);
        double* v     = velocity->getPointer(1);
        double* theta = dilatation->getPointer(0);
        double* omega = vorticity->getPointer(0);
        double* Omega = vorticity_magnitude->getPointer(0);
        std::vector<double*> F_node_x;
        std::vector<double*> F_node_y;
        F_node_x.reserve(d_num_eqn);
        F_node_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
            F_node_y.push_back(convective_flux_node[1]->getPointer(ei));
        }
        std::vector<double*> F_midpoint_x;
        std::vector<double*> F_midpoint_y;
        F_midpoint_x.reserve(d_num_eqn);
        F_midpoint_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
            F_midpoint_y.push_back(convective_flux_midpoint->getPointer(1, ei));
        }
        
        // Compute the magnitude of vorticity.
        for (int j = -num_conv_ghosts_1 + 1; j < interior_dim_1 + num_conv_ghosts_1 - 1; j++)
        {
            #ifdef __INTEL_COMPILER
            #pragma ivdep
            #endif
            for (int i = -num_conv_ghosts_0 + 1; i < interior_dim_0 + num_conv_ghosts_0 - 1; i++)
            {
                // Compute the linear indices.
                const int idx = (i + num_conv_ghosts_0) +
                    (j + num_conv_ghosts_1)*conv_ghostcell_dim_0;
                
                const int idx_vorticity = (i + num_subghosts_0_vorticity) +
                    (j + num_subghosts_1_vorticity)*subghostcell_dim_0_vorticity;
                
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
         * Declare temporary data containers for WENO interpolation.
         */
        
        std::vector<boost::shared_ptr<pdat::SideData<double> > > projection_variables;
        
        std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > > characteristic_variables;
        
        std::vector<boost::shared_ptr<pdat::SideData<double> > > characteristic_variables_minus;
        std::vector<boost::shared_ptr<pdat::SideData<double> > > characteristic_variables_plus;
        
        std::vector<boost::shared_ptr<pdat::SideData<double> > > primitive_variables_minus;
        std::vector<boost::shared_ptr<pdat::SideData<double> > > primitive_variables_plus;
        
        /*
         * Initialize temporary data containers for WENO interpolation.
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
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> W_array;
            W_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                W_array.push_back(characteristic_variables[m][ei]->getPointer(0));
            }
            
            double* W_L = characteristic_variables_minus[ei]->getPointer(0);
            double* W_R = characteristic_variables_plus[ei]->getPointer(0);
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                #ifdef __INTEL_COMPILER
                #pragma ivdep
                #endif
                for (int i = -1; i < interior_dim_0 + 2; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    #ifdef __INTEL_COMPILER
                    #pragma forceinline
                    #endif
                    performWENOInterpolation_minus(
                        W_L,
                        W_array,
                        idx_midpoint_x);
                }
            }
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                #ifdef __INTEL_COMPILER
                #pragma ivdep
                #endif
                for (int i = -1; i < interior_dim_0 + 2; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    #ifdef __INTEL_COMPILER
                    #pragma forceinline
                    #endif
                    performWENOInterpolation_plus(
                        W_R,
                        W_array,
                        idx_midpoint_x);
                }
            }
        }
        
        t_WENO_interpolation->stop();
        
        /*
         * Peform WENO interpolation in the y-direction.
         */
        
        t_WENO_interpolation->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> W_array;
            W_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                W_array.push_back(characteristic_variables[m][ei]->getPointer(1));
            }
            
            double* W_B = characteristic_variables_minus[ei]->getPointer(1);
            double* W_T = characteristic_variables_plus[ei]->getPointer(1);
            
            for (int j = -1; j < interior_dim_1 + 2; j++)
            {
                #ifdef __INTEL_COMPILER
                #pragma ivdep
                #endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    #ifdef __INTEL_COMPILER
                    #pragma forceinline
                    #endif
                    performWENOInterpolation_minus(
                        W_B,
                        W_array,
                        idx_midpoint_y);
                }
            }
            
            for (int j = -1; j < interior_dim_1 + 2; j++)
            {
                #ifdef __INTEL_COMPILER
                #pragma ivdep
                #endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    #ifdef __INTEL_COMPILER
                    #pragma forceinline
                    #endif
                    performWENOInterpolation_plus(
                        W_T,
                        W_array,
                        idx_midpoint_y);
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
        
        {
            std::vector<double*> V_L;
            std::vector<double*> V_R;
            V_L.reserve(d_num_eqn);
            V_R.reserve(d_num_eqn);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_L.push_back(primitive_variables_minus[ei]->getPointer(0));
                V_R.push_back(primitive_variables_plus[ei]->getPointer(0));
            }
            
            // Declare and initialize containers to store the references to the
            // WENO interpolated values.
            std::vector<boost::reference_wrapper<double> > V_L_ref;
            std::vector<boost::reference_wrapper<double> > V_R_ref;
            V_L_ref.reserve(d_num_eqn);
            V_R_ref.reserve(d_num_eqn);
            
            // Declare container to store the references to the mid-point flux.
            std::vector<boost::reference_wrapper<double> > F_midpoint_x_ref;
            F_midpoint_x_ref.reserve(d_num_eqn);
            
            // Declare container to store the references to the mid-point velocity.
            std::vector<boost::reference_wrapper<double> > vel_midpoint_x_ref;
            vel_midpoint_x_ref.reserve(d_dim.getValue());
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                for (int i = -1; i < interior_dim_0 + 2; i++)
                {
                    // Compute the linear index of the side.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                        
                    // Compute the average dilatation and magnitude of vorticity.
                    const int idx_L = (i - 1 + num_conv_ghosts_0) +
                        (j + num_conv_ghosts_1)*conv_ghostcell_dim_0;
                    
                    const int idx_R = (i + num_conv_ghosts_0) +
                        (j + num_conv_ghosts_1)*conv_ghostcell_dim_0;
                    
                    const int idx_dilatation_L = (i - 1 + num_subghosts_0_dilatation) +
                        (j + num_subghosts_1_dilatation)*subghostcell_dim_0_dilatation;
                    
                    const int idx_dilatation_R = (i + num_subghosts_0_dilatation) +
                        (j + num_subghosts_1_dilatation)*subghostcell_dim_0_dilatation;
                    
                    double theta_avg = 0.5*(theta[idx_dilatation_L] + theta[idx_dilatation_R]);
                    double Omega_avg = 0.5*(Omega[idx_L] + Omega[idx_R]);
                    
                    // Compute the Ducros-like shock sensor.
                    double s = -theta_avg/(fabs(theta_avg) + Omega_avg + EPSILON);
                    
                    // Initialize container that stores the references to the mid-point flux.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_midpoint_x_ref.push_back(boost::ref(F_midpoint_x[ei][idx_midpoint_x]));
                    }
                    
                    // Initialize container that stores the references to the mid-point velocity.
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        vel_midpoint_x_ref.push_back(
                            boost::ref(velocity_intercell->getPointer(0, di)[idx_midpoint_x]));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        V_L_ref.push_back(boost::ref(V_L[ei][idx_midpoint_x]));
                        V_R_ref.push_back(boost::ref(V_R[ei][idx_midpoint_x]));
                    }
                    
                    // Apply the Riemann solver.
                    if (s > 0.65)
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_midpoint_x_ref,
                                vel_midpoint_x_ref,
                                V_L_ref,
                                V_R_ref,
                                DIRECTION::X_DIRECTION,
                                RIEMANN_SOLVER::HLLC_HLL);
                    }
                    else
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_midpoint_x_ref,
                                vel_midpoint_x_ref,
                                V_L_ref,
                                V_R_ref,
                                DIRECTION::X_DIRECTION,
                                RIEMANN_SOLVER::HLLC);
                    }
                    
                    F_midpoint_x_ref.clear();
                    vel_midpoint_x_ref.clear();
                    V_L_ref.clear();
                    V_R_ref.clear();
                }
            }
        }
        
        t_Riemann_solver->stop();
        
        /*
         * Compute mid-point flux in the y-direction.
         */
        
        t_Riemann_solver->start();
        
        {
            std::vector<double*> V_B;
            std::vector<double*> V_T;
            V_B.reserve(d_num_eqn);
            V_T.reserve(d_num_eqn);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_B.push_back(primitive_variables_minus[ei]->getPointer(1));
                V_T.push_back(primitive_variables_plus[ei]->getPointer(1));
            }
            
            // Declare and initialize containers to store the references to the
            // WENO interpolated values.
            std::vector<boost::reference_wrapper<double> > V_B_ref;
            std::vector<boost::reference_wrapper<double> > V_T_ref;
            V_B_ref.reserve(d_num_eqn);
            V_T_ref.reserve(d_num_eqn);
            
            // Declare container to store the references to the mid-point flux.
            std::vector<boost::reference_wrapper<double> > F_midpoint_y_ref;
            F_midpoint_y_ref.reserve(d_num_eqn);
            
            // Declare container to store the references to the mid-point velocity.
            std::vector<boost::reference_wrapper<double> > vel_midpoint_y_ref;
            vel_midpoint_y_ref.reserve(d_dim.getValue());
            
            for (int j = -1; j < interior_dim_1 + 2; j++)
            {
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index of the side.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    // Compute the average dilatation and magnitude of vorticity.
                    const int idx_B = (i + num_conv_ghosts_0) +
                        (j - 1 + num_conv_ghosts_1)*conv_ghostcell_dim_0;
                    
                    const int idx_T = (i + num_conv_ghosts_0) +
                        (j + num_conv_ghosts_1)*conv_ghostcell_dim_0;
                    
                    const int idx_dilatation_B = (i + num_subghosts_0_dilatation) +
                        (j - 1 + num_subghosts_1_dilatation)*subghostcell_dim_0_dilatation;
                    
                    const int idx_dilatation_T = (i + num_subghosts_0_dilatation) +
                        (j + num_subghosts_1_dilatation)*subghostcell_dim_0_dilatation;
                    
                    double theta_avg = 0.5*(theta[idx_dilatation_B] + theta[idx_dilatation_T]);
                    double Omega_avg = 0.5*(Omega[idx_B] + Omega[idx_T]);
                    
                    // Compute the Ducros-like shock sensor.
                    double s = -theta_avg/(fabs(theta_avg) + Omega_avg + EPSILON);
                    
                    // Initialize container that stores the references to the mid-point flux.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_midpoint_y_ref.push_back(boost::ref(F_midpoint_y[ei][idx_midpoint_y]));
                    }
                    
                    // Initialize container that stores the references to the mid-point velocity.
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        vel_midpoint_y_ref.push_back(
                            boost::ref(velocity_intercell->getPointer(1, di)[idx_midpoint_y]));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        V_B_ref.push_back(boost::ref(V_B[ei][idx_midpoint_y]));
                        V_T_ref.push_back(boost::ref(V_T[ei][idx_midpoint_y]));
                    }
                    
                    // Apply the Riemann solver.
                    if (s > 0.65)
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_midpoint_y_ref,
                                vel_midpoint_y_ref,
                                V_B_ref,
                                V_T_ref,
                                DIRECTION::Y_DIRECTION,
                                RIEMANN_SOLVER::HLLC_HLL);
                    }
                    else
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_midpoint_y_ref,
                                vel_midpoint_y_ref,
                                V_B_ref,
                                V_T_ref,
                                DIRECTION::Y_DIRECTION,
                                RIEMANN_SOLVER::HLLC);
                    }
                    
                    F_midpoint_y_ref.clear();
                    vel_midpoint_y_ref.clear();
                    V_B_ref.clear();
                    V_T_ref.clear();
                }
            }
        }
        
        t_Riemann_solver->stop();
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        t_reconstruct_flux->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            double* F_face_x = convective_flux->getPointer(0, ei);
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                #ifdef __INTEL_COMPILER
                #pragma ivdep
                #endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    const int idx_midpoint_x_L = i +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    const int idx_midpoint_x_R = (i + 2) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    const int idx_node_L = (i - 1 + num_subghosts_0_convective_flux_x) +
                        (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                    
                    const int idx_node_R = (i + num_subghosts_0_convective_flux_x) +
                        (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                    
                    F_face_x[idx_face_x] = dt*(1.0/30.0*(F_midpoint_x[ei][idx_midpoint_x_R] +
                         F_midpoint_x[ei][idx_midpoint_x_L]) -
                        3.0/10.0*(F_node_x[ei][idx_node_R] +
                         F_node_x[ei][idx_node_L]) +
                        23.0/15.0*F_midpoint_x[ei][idx_midpoint_x]);
                }
            }
        }
        
        t_reconstruct_flux->stop();
        
        /*
         * Reconstruct the flux in the y-direction.
         */
        
        t_reconstruct_flux->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            double* F_face_y = convective_flux->getPointer(1, ei);
            
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
                #ifdef __INTEL_COMPILER
                #pragma ivdep
                #endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_y = j +
                        i*(interior_dim_1 + 1);
                    
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    const int idx_midpoint_y_B = (i + 1) +
                        j*(interior_dim_0 + 2);
                    
                    const int idx_midpoint_y_T = (i + 1) +
                        (j + 2)*(interior_dim_0 + 2);
                    
                    const int idx_node_B = (i + num_subghosts_0_convective_flux_y) +
                        (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                    
                    const int idx_node_T = (i + num_subghosts_0_convective_flux_y) +
                        (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                    
                    F_face_y[idx_face_y] = dt*(1.0/30.0*(F_midpoint_y[ei][idx_midpoint_y_T] +
                         F_midpoint_y[ei][idx_midpoint_y_B]) -
                        3.0/10.0*(F_node_y[ei][idx_node_T] +
                         F_node_y[ei][idx_node_B]) +
                        23.0/15.0*F_midpoint_y[ei][idx_midpoint_y]);
                }
            }
        }
        
        t_reconstruct_flux->stop();
        
        /*
         * Compute the source.
         */
        
        t_compute_source->start();
        
        const std::vector<EQN_FORM::TYPE> eqn_form = d_flow_model->getEquationsForm();
        
        double* u_midpoint_x = velocity_intercell->getPointer(0, 0);
        double* v_midpoint_y = velocity_intercell->getPointer(1, 1);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            if (eqn_form[ei] == EQN_FORM::ADVECTIVE)
            {
                double* S = source->getPointer(ei);
                
                const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
                const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
                
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #ifdef __INTEL_COMPILER
                    #pragma ivdep
                    #endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_velocity) +
                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                        
                        const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_velocity) +
                            (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                        
                        const int idx_cell_wghost_y_B = (i + num_subghosts_0_velocity) +
                            (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                        
                        const int idx_cell_wghost_y_T = (i + num_subghosts_0_velocity) +
                            (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity;
                        
                        const int idx_cell_nghost = i + j*interior_dim_0;
                        
                        const int idx_midpoint_x_LL = i + 
                            (j + 1)*(interior_dim_0 + 3);
                        
                        const int idx_midpoint_x_L = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3);
                        
                        const int idx_midpoint_x_R = (i + 2) +
                            (j + 1)*(interior_dim_0 + 3);
                        
                        const int idx_midpoint_x_RR = (i + 3) +
                            (j + 1)*(interior_dim_0 + 3);
                        
                        const int idx_midpoint_y_BB = (i + 1) +
                            j*(interior_dim_0 + 2);
                        
                        const int idx_midpoint_y_B = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2);
                        
                        const int idx_midpoint_y_T = (i + 1) +
                            (j + 2)*(interior_dim_0 + 2);
                        
                        const int idx_midpoint_y_TT = (i + 1) +
                            (j + 3)*(interior_dim_0 + 2);
                        
                        S[idx_cell_nghost] += dt*Q[ei][idx_cell_wghost]*(
                            (3.0/2.0*(u_midpoint_x[idx_midpoint_x_R] - u_midpoint_x[idx_midpoint_x_L]) -
                             3.0/10.0*(u[idx_cell_wghost_x_R] - u[idx_cell_wghost_x_L]) +
                             1.0/30.0*(u_midpoint_x[idx_midpoint_x_RR] - u_midpoint_x[idx_midpoint_x_LL]))/dx[0] +
                            (3.0/2.0*(v_midpoint_y[idx_midpoint_y_T] - v_midpoint_y[idx_midpoint_y_B]) -
                             3.0/10.0*(v[idx_cell_wghost_y_T] - v[idx_cell_wghost_y_B]) +
                             1.0/30.0*(v_midpoint_y[idx_midpoint_y_TT] - v_midpoint_y[idx_midpoint_y_BB]))/dx[1]);
                    }
                }
            }
        }
        
        t_compute_source->stop();
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_conv_ghosts_0 = d_num_conv_ghosts[0];
        const int num_conv_ghosts_1 = d_num_conv_ghosts[1];
        const int num_conv_ghosts_2 = d_num_conv_ghosts[2];
        
        const int conv_ghostcell_dim_0 = conv_ghostcell_dims[0];
        const int conv_ghostcell_dim_1 = conv_ghostcell_dims[1];
        
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
        
        const int num_subghosts_0_velocity = num_subghosts_velocity[0];
        const int num_subghosts_1_velocity = num_subghosts_velocity[1];
        const int num_subghosts_2_velocity = num_subghosts_velocity[2];
        const int subghostcell_dim_0_velocity = subghostcell_dims_velocity[0];
        const int subghostcell_dim_1_velocity = subghostcell_dims_velocity[1];
        
        const int num_subghosts_0_dilatation = num_subghosts_dilatation[0];
        const int num_subghosts_1_dilatation = num_subghosts_dilatation[1];
        const int num_subghosts_2_dilatation = num_subghosts_dilatation[2];
        const int subghostcell_dim_0_dilatation = subghostcell_dims_dilatation[0];
        const int subghostcell_dim_1_dilatation = subghostcell_dims_dilatation[1];
        
        const int num_subghosts_0_vorticity = num_subghosts_vorticity[0];
        const int num_subghosts_1_vorticity = num_subghosts_vorticity[1];
        const int num_subghosts_2_vorticity = num_subghosts_vorticity[2];
        const int subghostcell_dim_0_vorticity = subghostcell_dims_vorticity[0];
        const int subghostcell_dim_1_vorticity = subghostcell_dims_vorticity[1];
        
        const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
        const int num_subghosts_1_convective_flux_x = num_subghosts_convective_flux_x[1];
        const int num_subghosts_2_convective_flux_x = num_subghosts_convective_flux_x[2];
        const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
        const int subghostcell_dim_1_convective_flux_x = subghostcell_dims_convective_flux_x[1];
        
        const int num_subghosts_0_convective_flux_y = num_subghosts_convective_flux_y[0];
        const int num_subghosts_1_convective_flux_y = num_subghosts_convective_flux_y[1];
        const int num_subghosts_2_convective_flux_y = num_subghosts_convective_flux_y[2];
        const int subghostcell_dim_0_convective_flux_y = subghostcell_dims_convective_flux_y[0];
        const int subghostcell_dim_1_convective_flux_y = subghostcell_dims_convective_flux_y[1];
        
        const int num_subghosts_0_convective_flux_z = num_subghosts_convective_flux_z[0];
        const int num_subghosts_1_convective_flux_z = num_subghosts_convective_flux_z[1];
        const int num_subghosts_2_convective_flux_z = num_subghosts_convective_flux_z[2];
        const int subghostcell_dim_0_convective_flux_z = subghostcell_dims_convective_flux_z[0];
        const int subghostcell_dim_1_convective_flux_z = subghostcell_dims_convective_flux_z[1];
        
        double* u     = velocity->getPointer(0);
        double* v     = velocity->getPointer(1);
        double* w     = velocity->getPointer(2);
        double* theta = dilatation->getPointer(0);
        double* omega_x = vorticity->getPointer(0);
        double* omega_y = vorticity->getPointer(1);
        double* omega_z = vorticity->getPointer(2);
        double* Omega = vorticity_magnitude->getPointer(0);
        std::vector<double*> F_node_x;
        std::vector<double*> F_node_y;
        std::vector<double*> F_node_z;
        F_node_x.reserve(d_num_eqn);
        F_node_y.reserve(d_num_eqn);
        F_node_z.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
            F_node_y.push_back(convective_flux_node[1]->getPointer(ei));
            F_node_z.push_back(convective_flux_node[2]->getPointer(ei));
        }
        std::vector<double*> F_midpoint_x;
        std::vector<double*> F_midpoint_y;
        std::vector<double*> F_midpoint_z;
        F_midpoint_x.reserve(d_num_eqn);
        F_midpoint_y.reserve(d_num_eqn);
        F_midpoint_z.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
            F_midpoint_y.push_back(convective_flux_midpoint->getPointer(1, ei));
            F_midpoint_z.push_back(convective_flux_midpoint->getPointer(2, ei));
        }
        
        // Compute the magnitude of vorticity.
        for (int k = -num_conv_ghosts_2 + 1; k < interior_dim_2 + num_conv_ghosts_2 - 1; k++)
        {
            for (int j = -num_conv_ghosts_1 + 1; j < interior_dim_1 + num_conv_ghosts_1 - 1; j++)
            {
                #ifdef __INTEL_COMPILER
                #pragma ivdep
                #endif
                for (int i = -num_conv_ghosts_0 + 1; i < interior_dim_0 + num_conv_ghosts_0 - 1; i++)
                {
                    // Compute the linear indices
                    const int idx = (i + num_conv_ghosts_0) +
                        (j + num_conv_ghosts_1)*conv_ghostcell_dim_0 +
                        (k + num_conv_ghosts_2)*conv_ghostcell_dim_0*
                            conv_ghostcell_dim_1;
                    
                    const int idx_vorticity = (i + num_subghosts_0_vorticity) +
                        (j + num_subghosts_1_vorticity)*subghostcell_dim_0_vorticity +
                        (k + num_subghosts_2_vorticity)*subghostcell_dim_0_vorticity*
                            subghostcell_dim_1_vorticity;
                    
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
         * Declare temporary data containers for WENO interpolation.
         */
        
        std::vector<boost::shared_ptr<pdat::SideData<double> > > projection_variables;
        
        std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > > characteristic_variables;
        
        std::vector<boost::shared_ptr<pdat::SideData<double> > > characteristic_variables_minus;
        std::vector<boost::shared_ptr<pdat::SideData<double> > > characteristic_variables_plus;
        
        std::vector<boost::shared_ptr<pdat::SideData<double> > > primitive_variables_minus;
        std::vector<boost::shared_ptr<pdat::SideData<double> > > primitive_variables_plus;
        
        /*
         * Initialize temporary data containers for WENO interpolation.
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
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> W_array;
            W_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                W_array.push_back(characteristic_variables[m][ei]->getPointer(0));
            }
            
            double* W_L = characteristic_variables_minus[ei]->getPointer(0);
            double* W_R = characteristic_variables_plus[ei]->getPointer(0);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #ifdef __INTEL_COMPILER
                    #pragma ivdep
                    #endif
                    for (int i = -1; i < interior_dim_0 + 2; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        #ifdef __INTEL_COMPILER
                        #pragma forceinline
                        #endif
                        performWENOInterpolation_minus(
                            W_L,
                            W_array,
                            idx_midpoint_x);
                    }
                }
            }
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #ifdef __INTEL_COMPILER
                    #pragma ivdep
                    #endif
                    for (int i = -1; i < interior_dim_0 + 2; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        #ifdef __INTEL_COMPILER
                        #pragma forceinline
                        #endif
                        performWENOInterpolation_plus(
                            W_R,
                            W_array,
                            idx_midpoint_x);
                    }
                }
            }
        }
        
        t_WENO_interpolation->stop();
        
        /*
         * Peform WENO interpolation in the y-direction.
         */
        
        t_WENO_interpolation->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> W_array;
            W_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                W_array.push_back(characteristic_variables[m][ei]->getPointer(1));
            }
            
            double* W_B = characteristic_variables_minus[ei]->getPointer(1);
            double* W_T = characteristic_variables_plus[ei]->getPointer(1);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -1; j < interior_dim_1 + 2; j++)
                {
                    #ifdef __INTEL_COMPILER
                    #pragma ivdep
                    #endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        #ifdef __INTEL_COMPILER
                        #pragma forceinline
                        #endif
                        performWENOInterpolation_minus(
                            W_B,
                            W_array,
                            idx_midpoint_y);
                    }
                }
            }
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -1; j < interior_dim_1 + 2; j++)
                {
                    #ifdef __INTEL_COMPILER
                    #pragma ivdep
                    #endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        #ifdef __INTEL_COMPILER
                        #pragma forceinline
                        #endif
                        performWENOInterpolation_plus(
                            W_T,
                            W_array,
                            idx_midpoint_y);
                    }
                }
            }
        }
        
        t_WENO_interpolation->stop();
        
        /*
         * Peform WENO interpolation in the z-direction.
         */
        
        t_WENO_interpolation->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> W_array;
            W_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                W_array.push_back(characteristic_variables[m][ei]->getPointer(2));
            }
            
            double* W_B = characteristic_variables_minus[ei]->getPointer(2);
            double* W_F = characteristic_variables_plus[ei]->getPointer(2);
            
            for (int k = -1; k < interior_dim_2 + 2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #ifdef __INTEL_COMPILER
                    #pragma ivdep
                    #endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        #ifdef __INTEL_COMPILER
                        #pragma forceinline
                        #endif
                        performWENOInterpolation_minus(
                            W_B,
                            W_array,
                            idx_midpoint_z);
                    }
                }
            }
            
            for (int k = -1; k < interior_dim_2 + 2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #ifdef __INTEL_COMPILER
                    #pragma ivdep
                    #endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        #ifdef __INTEL_COMPILER
                        #pragma forceinline
                        #endif
                        performWENOInterpolation_plus(
                            W_F,
                            W_array,
                            idx_midpoint_z);
                    }
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
        
        {
            std::vector<double*> V_L;
            std::vector<double*> V_R;
            V_L.reserve(d_num_eqn);
            V_R.reserve(d_num_eqn);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_L.push_back(primitive_variables_minus[ei]->getPointer(0));
                V_R.push_back(primitive_variables_plus[ei]->getPointer(0));
            }
            
            // Declare and initialize containers to store the references to the
            // WENO interpolated values.
            std::vector<boost::reference_wrapper<double> > V_L_ref;
            std::vector<boost::reference_wrapper<double> > V_R_ref;
            V_L_ref.reserve(d_num_eqn);
            V_R_ref.reserve(d_num_eqn);
            
            // Declare container to store the references to the mid-point flux.
            std::vector<boost::reference_wrapper<double> > F_midpoint_x_ref;
            F_midpoint_x_ref.reserve(d_num_eqn);
            
            // Declare container to store the references to the mid-point velocity.
            std::vector<boost::reference_wrapper<double> > vel_midpoint_x_ref;
            vel_midpoint_x_ref.reserve(d_dim.getValue());
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    for (int i = -1; i < interior_dim_0 + 2; i++)
                    {
                        // Compute the linear index of the side.
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        // Compute the average dilatation and magnitude of vorticity.
                        const int idx_L = (i - 1 + num_conv_ghosts_0) +
                            (j + num_conv_ghosts_1)*conv_ghostcell_dim_0 +
                            (k + num_conv_ghosts_2)*conv_ghostcell_dim_0*
                                conv_ghostcell_dim_1;
                            
                        
                        const int idx_R = (i + num_conv_ghosts_0) +
                            (j + num_conv_ghosts_1)*conv_ghostcell_dim_0 +
                            (k + num_conv_ghosts_2)*conv_ghostcell_dim_0*
                                conv_ghostcell_dim_1;
                        
                        const int idx_dilatation_L = (i - 1 + num_subghosts_0_dilatation) +
                            (j + num_subghosts_1_dilatation)*subghostcell_dim_0_dilatation +
                            (k + num_subghosts_2_dilatation)*subghostcell_dim_0_dilatation*
                                subghostcell_dim_1_dilatation;
                        
                        const int idx_dilatation_R = (i + num_subghosts_0_dilatation) +
                            (j + num_subghosts_1_dilatation)*subghostcell_dim_0_dilatation +
                            (k + num_subghosts_2_dilatation)*subghostcell_dim_0_dilatation*
                                subghostcell_dim_1_dilatation;
                        
                        double theta_avg = 0.5*(theta[idx_dilatation_L] + theta[idx_dilatation_R]);
                        double Omega_avg = 0.5*(Omega[idx_L] + Omega[idx_R]);
                        
                        // Compute the Ducros-like shock sensor.
                        double s = -theta_avg/(fabs(theta_avg) + Omega_avg + EPSILON);
                        
                        // Initialize container that stores the references to the mid-point flux.
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            F_midpoint_x_ref.push_back(boost::ref(F_midpoint_x[ei][idx_midpoint_x]));
                        }
                        
                        // Initialize container that stores the references to the mid-point velocity.
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            vel_midpoint_x_ref.push_back(
                                boost::ref(velocity_intercell->getPointer(0, di)[idx_midpoint_x]));
                        }
                        
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            V_L_ref.push_back(boost::ref(V_L[ei][idx_midpoint_x]));
                            V_R_ref.push_back(boost::ref(V_R[ei][idx_midpoint_x]));
                        }
                        
                        // Apply the Riemann solver.
                        if (s > 0.65)
                        {
                            d_flow_model->
                                computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                    F_midpoint_x_ref,
                                    vel_midpoint_x_ref,
                                    V_L_ref,
                                    V_R_ref,
                                    DIRECTION::X_DIRECTION,
                                    RIEMANN_SOLVER::HLLC_HLL);
                        }
                        else
                        {
                            d_flow_model->
                                computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                    F_midpoint_x_ref,
                                    vel_midpoint_x_ref,
                                    V_L_ref,
                                    V_R_ref,
                                    DIRECTION::X_DIRECTION,
                                    RIEMANN_SOLVER::HLLC);
                        }
                        
                        F_midpoint_x_ref.clear();
                        vel_midpoint_x_ref.clear();
                        V_L_ref.clear();
                        V_R_ref.clear();
                    }
                }
            }
        }
        
        t_Riemann_solver->stop();
        
        /*
         * Compute mid-point flux in the y-direction.
         */
        
        t_Riemann_solver->start();
        
        {
            std::vector<double*> V_B;
            std::vector<double*> V_T;
            V_B.reserve(d_num_eqn);
            V_T.reserve(d_num_eqn);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_B.push_back(primitive_variables_minus[ei]->getPointer(1));
                V_T.push_back(primitive_variables_plus[ei]->getPointer(1));
            }
            
            // Declare and initialize containers to store the references to the
            // WENO interpolated values.
            std::vector<boost::reference_wrapper<double> > V_B_ref;
            std::vector<boost::reference_wrapper<double> > V_T_ref;
            V_B_ref.reserve(d_num_eqn);
            V_T_ref.reserve(d_num_eqn);
            
            // Declare container to store the references to the mid-point flux.
            std::vector<boost::reference_wrapper<double> > F_midpoint_y_ref;
            F_midpoint_y_ref.reserve(d_num_eqn);
            
            // Declare container to store the references to the mid-point velocity.
            std::vector<boost::reference_wrapper<double> > vel_midpoint_y_ref;
            vel_midpoint_y_ref.reserve(d_dim.getValue());
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -1; j < interior_dim_1 + 2; j++)
                {
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the side.
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        // Compute the average dilatation and magnitude of vorticity.
                        const int idx_B = (i + num_conv_ghosts_0) +
                            (j - 1 + num_conv_ghosts_1)*conv_ghostcell_dim_0 +
                            (k + num_conv_ghosts_2)*conv_ghostcell_dim_0*
                                conv_ghostcell_dim_1;
                        
                        const int idx_T = (i + num_conv_ghosts_0) +
                            (j + num_conv_ghosts_1)*conv_ghostcell_dim_0 +
                            (k + num_conv_ghosts_2)*conv_ghostcell_dim_0*
                                conv_ghostcell_dim_1;
                        
                        const int idx_dilatation_B = (i + num_subghosts_0_dilatation) +
                            (j - 1 + num_subghosts_1_dilatation)*subghostcell_dim_0_dilatation +
                            (k + num_subghosts_2_dilatation)*subghostcell_dim_0_dilatation*
                                subghostcell_dim_1_dilatation;
                        
                        const int idx_dilatation_T = (i + num_subghosts_0_dilatation) +
                            (j + num_subghosts_1_dilatation)*subghostcell_dim_0_dilatation +
                            (k + num_subghosts_2_dilatation)*subghostcell_dim_0_dilatation*
                                subghostcell_dim_1_dilatation;
                        
                        double theta_avg = 0.5*(theta[idx_dilatation_B] + theta[idx_dilatation_T]);
                        double Omega_avg = 0.5*(Omega[idx_B] + Omega[idx_T]);
                        
                        // Compute the Ducros-like shock sensor.
                        double s = -theta_avg/(fabs(theta_avg) + Omega_avg + EPSILON);
                        
                        // Initialize container that stores the references to the mid-point flux.
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            F_midpoint_y_ref.push_back(boost::ref(F_midpoint_y[ei][idx_midpoint_y]));
                        }
                        
                        // Initialize container that stores the references to the mid-point velocity.
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            vel_midpoint_y_ref.push_back(
                                boost::ref(velocity_intercell->getPointer(1, di)[idx_midpoint_y]));
                        }
                        
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            V_B_ref.push_back(boost::ref(V_B[ei][idx_midpoint_y]));
                            V_T_ref.push_back(boost::ref(V_T[ei][idx_midpoint_y]));
                        }
                        
                        // Apply the Riemann solver.
                        if (s > 0.65)
                        {
                            d_flow_model->
                                computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                    F_midpoint_y_ref,
                                    vel_midpoint_y_ref,
                                    V_B_ref,
                                    V_T_ref,
                                    DIRECTION::Y_DIRECTION,
                                    RIEMANN_SOLVER::HLLC_HLL);
                        }
                        else
                        {
                            d_flow_model->
                                computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                    F_midpoint_y_ref,
                                    vel_midpoint_y_ref,
                                    V_B_ref,
                                    V_T_ref,
                                    DIRECTION::Y_DIRECTION,
                                    RIEMANN_SOLVER::HLLC);
                        }
                        
                        F_midpoint_y_ref.clear();
                        vel_midpoint_y_ref.clear();
                        V_B_ref.clear();
                        V_T_ref.clear();
                    }
                }
            }
        }
        
        t_Riemann_solver->stop();
        
        /*
         * Compute mid-point flux in the z-direction.
         */
        
        t_Riemann_solver->start();
        
        {
            std::vector<double*> V_B;
            std::vector<double*> V_F;
            V_B.reserve(d_num_eqn);
            V_F.reserve(d_num_eqn);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_B.push_back(primitive_variables_minus[ei]->getPointer(2));
                V_F.push_back(primitive_variables_plus[ei]->getPointer(2));
            }
            
            // Declare and initialize containers to store the references to the
            // WENO interpolated values.
            std::vector<boost::reference_wrapper<double> > V_B_ref;
            std::vector<boost::reference_wrapper<double> > V_F_ref;
            V_B_ref.reserve(d_num_eqn);
            V_F_ref.reserve(d_num_eqn);
            
            // Declare container to store the references to the mid-point flux.
            std::vector<boost::reference_wrapper<double> > F_midpoint_z_ref;
            F_midpoint_z_ref.reserve(d_num_eqn);
            
            // Declare container to store the references to the mid-point velocity.
            std::vector<boost::reference_wrapper<double> > vel_midpoint_z_ref;
            vel_midpoint_z_ref.reserve(d_dim.getValue());
            
            for (int k = -1; k < interior_dim_2 + 2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the side.
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        // Compute the average dilatation and magnitude of vorticity.
                        const int idx_B = (i + num_conv_ghosts_0) +
                            (j + num_conv_ghosts_1)*conv_ghostcell_dim_0 +
                            (k - 1 + num_conv_ghosts_2)*conv_ghostcell_dim_0*
                                conv_ghostcell_dim_1;
                        
                        const int idx_F = (i + num_conv_ghosts_0) +
                            (j + num_conv_ghosts_1)*conv_ghostcell_dim_0 +
                            (k + num_conv_ghosts_2)*conv_ghostcell_dim_0*
                                conv_ghostcell_dim_1;
                        
                        const int idx_dilatation_B = (i + num_subghosts_0_dilatation) +
                            (j + num_subghosts_1_dilatation)*subghostcell_dim_0_dilatation +
                            (k - 1 + num_subghosts_2_dilatation)*subghostcell_dim_0_dilatation*
                                subghostcell_dim_1_dilatation;
                        
                        const int idx_dilatation_F = (i + num_subghosts_0_dilatation) +
                            (j + num_subghosts_1_dilatation)*subghostcell_dim_0_dilatation +
                            (k + num_subghosts_2_dilatation)*subghostcell_dim_0_dilatation*
                                subghostcell_dim_1_dilatation;
                        
                        double theta_avg = 0.5*(theta[idx_dilatation_B] + theta[idx_dilatation_F]);
                        double Omega_avg = 0.5*(Omega[idx_B] + Omega[idx_F]);
                        
                        // Compute the Ducros-like shock sensor.
                        double s = -theta_avg/(fabs(theta_avg) + Omega_avg + EPSILON);
                        
                        // Initialize container that stores the references to the mid-point flux.
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            F_midpoint_z_ref.push_back(boost::ref(F_midpoint_z[ei][idx_midpoint_z]));
                        }
                        
                        // Initialize container that stores the references to the mid-point velocity.
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            vel_midpoint_z_ref.push_back(
                                boost::ref(velocity_intercell->getPointer(2, di)[idx_midpoint_z]));
                        }
                        
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            V_B_ref.push_back(boost::ref(V_B[ei][idx_midpoint_z]));
                            V_F_ref.push_back(boost::ref(V_F[ei][idx_midpoint_z]));
                        }
                        
                        // Apply the Riemann solver.
                        if (s > 0.65)
                        {
                            d_flow_model->
                                computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                    F_midpoint_z_ref,
                                    vel_midpoint_z_ref,
                                    V_B_ref,
                                    V_F_ref,
                                    DIRECTION::Z_DIRECTION,
                                    RIEMANN_SOLVER::HLLC_HLL);
                        }
                        else
                        {
                            d_flow_model->
                                computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                    F_midpoint_z_ref,
                                    vel_midpoint_z_ref,
                                    V_B_ref,
                                    V_F_ref,
                                    DIRECTION::Z_DIRECTION,
                                    RIEMANN_SOLVER::HLLC);
                        }
                        
                        F_midpoint_z_ref.clear();
                        vel_midpoint_z_ref.clear();
                        V_B_ref.clear();
                        V_F_ref.clear();
                    }
                }
            }
        }
        
        t_Riemann_solver->stop();
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        t_reconstruct_flux->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            double* F_face_x = convective_flux->getPointer(0, ei);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #ifdef __INTEL_COMPILER
                    #pragma ivdep
                    #endif
                    for (int i = 0; i < interior_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dim_0 + 1) +
                            k*(interior_dim_0 + 1)*
                                interior_dim_1;
                        
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        const int idx_midpoint_x_L = i +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        const int idx_midpoint_x_R = (i + 2) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        const int idx_node_L = (i - 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                            (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                subghostcell_dim_1_convective_flux_x;
                        
                        const int idx_node_R = (i + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                            (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                subghostcell_dim_1_convective_flux_x;
                        
                        F_face_x[idx_face_x] = dt*(1.0/30.0*(F_midpoint_x[ei][idx_midpoint_x_R] +
                             F_midpoint_x[ei][idx_midpoint_x_L]) -
                            3.0/10.0*(F_node_x[ei][idx_node_R] +
                             F_node_x[ei][idx_node_L]) +
                            23.0/15.0*F_midpoint_x[ei][idx_midpoint_x]);
                    }
                }
            }
        }
        
        t_reconstruct_flux->stop();
        
        /*
         * Reconstruct the flux in the y-direction.
         */
        
        t_reconstruct_flux->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            double* F_face_y = convective_flux->getPointer(1, ei);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
                    #ifdef __INTEL_COMPILER
                    #pragma ivdep
                    #endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = j +
                            k*(interior_dim_1 + 1) +
                            i*(interior_dim_1 + 1)*
                                interior_dim_2;
                        
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        const int idx_midpoint_y_B = (i + 1) +
                            j*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        const int idx_midpoint_y_T = (i + 1) +
                            (j + 2)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        const int idx_node_B = (i + num_subghosts_0_convective_flux_y) +
                            (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                            (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                subghostcell_dim_1_convective_flux_y;
                        
                        const int idx_node_T = (i + num_subghosts_0_convective_flux_y) +
                            (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                            (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                subghostcell_dim_1_convective_flux_y;
                        
                        F_face_y[idx_face_y] = dt*(1.0/30.0*(F_midpoint_y[ei][idx_midpoint_y_T] +
                             F_midpoint_y[ei][idx_midpoint_y_B]) -
                            3.0/10.0*(F_node_y[ei][idx_node_T] +
                             F_node_y[ei][idx_node_B]) +
                            23.0/15.0*F_midpoint_y[ei][idx_midpoint_y]);
                    }
                }
            }
        }
        
        t_reconstruct_flux->stop();
        
        /*
         * Reconstruct the flux in the z-direction.
         */
        
        t_reconstruct_flux->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            double* F_face_z = convective_flux->getPointer(2, ei);
            
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #ifdef __INTEL_COMPILER
                    #pragma ivdep
                    #endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_z = k +
                            i*(interior_dims[2] + 1) +
                            j*(interior_dims[2] + 1)*interior_dims[0];
                        
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const int idx_midpoint_z_B = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            k*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const int idx_midpoint_z_F = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 2)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const int idx_node_B = (i + num_subghosts_0_convective_flux_z) +
                            (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                            (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                subghostcell_dim_1_convective_flux_z;
                        
                        const int idx_node_F = (i + num_subghosts_0_convective_flux_z) +
                            (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                            (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                subghostcell_dim_1_convective_flux_z;
                        
                        F_face_z[idx_face_z] = dt*(1.0/30.0*(F_midpoint_z[ei][idx_midpoint_z_F] +
                             F_midpoint_z[ei][idx_midpoint_z_B]) -
                            3.0/10.0*(F_node_z[ei][idx_node_F] +
                             F_node_z[ei][idx_node_B]) +
                            23.0/15.0*F_midpoint_z[ei][idx_midpoint_z]);
                    }
                }
            }
        }
        
        t_reconstruct_flux->stop();
        
        /*
         * Compute the source.
         */
        
        t_compute_source->start();
        
        const std::vector<EQN_FORM::TYPE> eqn_form = d_flow_model->getEquationsForm();
        
        double* u_midpoint_x = velocity_intercell->getPointer(0, 0);
        double* v_midpoint_y = velocity_intercell->getPointer(1, 1);
        double* w_midpoint_z = velocity_intercell->getPointer(2, 2);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            if (eqn_form[ei] == EQN_FORM::ADVECTIVE)
            {
                double* S = source->getPointer(ei);
                
                const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
                const int num_subghosts_2_conservative_var = num_subghosts_conservative_var[ei][2];
                const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
                const int subghostcell_dim_1_conservative_var = subghostcell_dims_conservative_var[ei][1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        #ifdef __INTEL_COMPILER
                        #pragma ivdep
                        #endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices. 
                            const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                    subghostcell_dim_1_velocity;
                            
                            const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                    subghostcell_dim_1_velocity;
                            
                            const int idx_cell_wghost_y_B = (i + num_subghosts_0_velocity) +
                                (j - 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                    subghostcell_dim_1_velocity;
                            
                            const int idx_cell_wghost_y_T = (i + num_subghosts_0_velocity) +
                                (j + 1 + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                (k + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                    subghostcell_dim_1_velocity;
                            
                            const int idx_cell_wghost_z_B = (i + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                (k - 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                    subghostcell_dim_1_velocity;
                            
                            const int idx_cell_wghost_z_F = (i + num_subghosts_0_velocity) +
                                (j + num_subghosts_1_velocity)*subghostcell_dim_0_velocity +
                                (k + 1 + num_subghosts_2_velocity)*subghostcell_dim_0_velocity*
                                    subghostcell_dim_1_velocity;
                            
                            const int idx_cell_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_midpoint_x_LL = i +
                                (j + 1)*(interior_dim_0 + 3) +
                                (k + 1)*(interior_dim_0 + 3)*
                                    (interior_dim_1 + 2);
                            
                            const int idx_midpoint_x_L = (i + 1) +
                                (j + 1)*(interior_dim_0 + 3) +
                                (k + 1)*(interior_dim_0 + 3)*
                                    (interior_dim_1 + 2);
                            
                            const int idx_midpoint_x_R = (i + 2) +
                                (j + 1)*(interior_dim_0 + 3) +
                                (k + 1)*(interior_dim_0 + 3)*
                                    (interior_dim_1 + 2);
                            
                            const int idx_midpoint_x_RR = (i + 3) +
                                (j + 1)*(interior_dim_0 + 3) +
                                (k + 1)*(interior_dim_0 + 3)*
                                    (interior_dim_1 + 2);
                            
                            const int idx_midpoint_y_BB = (i + 1) +
                                j*(interior_dim_0 + 2) +
                                (k + 1)*(interior_dim_0 + 2)*
                                    (interior_dim_1 + 3);
                            
                            const int idx_midpoint_y_B = (i + 1) +
                                (j + 1)*(interior_dim_0 + 2) +
                                (k + 1)*(interior_dim_0 + 2)*
                                    (interior_dim_1 + 3);
                            
                            const int idx_midpoint_y_T = (i + 1) +
                                (j + 2)*(interior_dim_0 + 2) +
                                (k + 1)*(interior_dim_0 + 2)*
                                    (interior_dim_1 + 3);
                            
                            const int idx_midpoint_y_TT = (i + 1) +
                                (j + 3)*(interior_dim_0 + 2) +
                                (k + 1)*(interior_dim_0 + 2)*
                                    (interior_dim_1 + 3);
                            
                            const int idx_midpoint_z_BB = (i + 1) +
                                (j + 1)*(interior_dim_0 + 2) +
                                k*(interior_dim_0 + 2)*
                                    (interior_dim_1 + 2);
                            
                            const int idx_midpoint_z_B = (i + 1) +
                                (j + 1)*(interior_dim_0 + 2) +
                                (k + 1)*(interior_dim_0 + 2)*
                                    (interior_dim_1 + 2);
                            
                            const int idx_midpoint_z_F = (i + 1) +
                                (j + 1)*(interior_dim_0 + 2) +
                                (k + 2)*(interior_dim_0 + 2)*
                                    (interior_dim_1 + 2);
                            
                            const int idx_midpoint_z_FF = (i + 1) +
                                (j + 1)*(interior_dim_0 + 2) +
                                (k + 3)*(interior_dim_0 + 2)*
                                    (interior_dim_1 + 2);
                            
                            S[idx_cell_nghost] += dt*Q[ei][idx_cell_wghost]*(
                                (3.0/2.0*(u_midpoint_x[idx_midpoint_x_R] - u_midpoint_x[idx_midpoint_x_L]) -
                                 3.0/10.0*(u[idx_cell_wghost_x_R] - u[idx_cell_wghost_x_L]) +
                                 1.0/30.0*(u_midpoint_x[idx_midpoint_x_RR] - u_midpoint_x[idx_midpoint_x_LL]))/dx[0] +
                                (3.0/2.0*(v_midpoint_y[idx_midpoint_y_T] - v_midpoint_y[idx_midpoint_y_B]) -
                                 3.0/10.0*(v[idx_cell_wghost_y_T] - v[idx_cell_wghost_y_B]) +
                                 1.0/30.0*(v_midpoint_y[idx_midpoint_y_TT] - v_midpoint_y[idx_midpoint_y_BB]))/dx[1] +
                                (3.0/2.0*(w_midpoint_z[idx_midpoint_z_F] - w_midpoint_z[idx_midpoint_z_B]) -
                                 3.0/10.0*(w[idx_cell_wghost_z_F] - w[idx_cell_wghost_z_B]) +
                                 1.0/30.0*(w_midpoint_z[idx_midpoint_z_FF] - w_midpoint_z[idx_midpoint_z_BB]))/dx[2]);
                        }
                    }
                }
            }
        }
        
        t_compute_source->stop();
        
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
    const double& alpha_tau = d_constant_alpha_tau;
    
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
        
        if(fabs(tau_6/(beta_avg + EPSILON)) > alpha_tau)
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
        
        if (fabs(tau_6_tilde/(beta_tilde_avg + EPSILON)) > alpha_tau)
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
    double* beta_0,
    double* beta_1,
    double* beta_2,
    double* beta_3,
    const std::vector<double*>& U_array,
    const int& idx_side)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.size()) == 6);
#endif
    
    *beta_0 = 1.0/3.0*(U_array[0][idx_side]*(4.0*U_array[0][idx_side] - 19.0*U_array[1][idx_side] +
         11.0*U_array[2][idx_side]) + U_array[1][idx_side]*(25.0*U_array[1][idx_side] -
         31.0*U_array[2][idx_side]) + 10.0*U_array[2][idx_side]*U_array[2][idx_side]);
    
    *beta_1 = 1.0/3.0*(U_array[1][idx_side]*(4.0*U_array[1][idx_side] - 13.0*U_array[2][idx_side] +
         5.0*U_array[3][idx_side]) + 13.0*U_array[2][idx_side]*(U_array[2][idx_side] -
         U_array[3][idx_side]) + 4.0*U_array[3][idx_side]*U_array[3][idx_side]);
    
    *beta_2 = 1.0/3.0*(U_array[2][idx_side]*(10.0*U_array[2][idx_side] - 31.0*U_array[3][idx_side] +
         11.0*U_array[4][idx_side]) + U_array[3][idx_side]*(25.0*U_array[3][idx_side] -
         19.0*U_array[4][idx_side]) + 4.0*U_array[4][idx_side]*U_array[4][idx_side]);
    
    *beta_3 = 1.0/232243200.0*(U_array[0][idx_side]*(525910327.0*U_array[0][idx_side] -
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
    double* beta_tilde_0,
    double* beta_tilde_1,
    double* beta_tilde_2,
    double* beta_tilde_3,
    const std::vector<double*>& U_array,
    const int& idx_side)
{
#ifdef DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(U_array.size()) == 6);
#endif
    
    *beta_tilde_0 = 1.0/3.0*(U_array[5][idx_side]*(4.0*U_array[5][idx_side] - 19.0*U_array[4][idx_side] +
         11.0*U_array[3][idx_side]) + U_array[4][idx_side]*(25.0*U_array[4][idx_side] -
         31.0*U_array[3][idx_side]) + 10.0*U_array[3][idx_side]*U_array[3][idx_side]);
    
    *beta_tilde_1 = 1.0/3.0*(U_array[4][idx_side]*(4.0*U_array[4][idx_side] - 13.0*U_array[3][idx_side] +
         5.0*U_array[2][idx_side]) + 13.0*U_array[3][idx_side]*(U_array[3][idx_side] -
         U_array[2][idx_side]) + 4.0*U_array[2][idx_side]*U_array[2][idx_side]);
    
    *beta_tilde_2 = 1.0/3.0*(U_array[3][idx_side]*(10.0*U_array[3][idx_side] - 31.0*U_array[2][idx_side] +
         11.0*U_array[1][idx_side]) + U_array[2][idx_side]*(25.0*U_array[2][idx_side] -
         19.0*U_array[1][idx_side]) + 4.0*U_array[1][idx_side]*U_array[1][idx_side]);
    
    *beta_tilde_3 = 1.0/232243200.0*(U_array[5][idx_side]*(525910327.0*U_array[5][idx_side] -
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
ConvectiveFluxReconstructorWCNS6_Test::performWENOInterpolation_minus(
    double* U_minus,
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
    const double& alpha_tau  = d_constant_alpha_tau;
    
    /*
     * Compute sigma.
     */
    
    double sigma;
    
    #ifdef __INTEL_COMPILER
    #pragma forceinline
    #endif
    computeSigma(&sigma, U_array, idx_side);
    
    /*
     * Compute beta's.
     */
    
    double beta_0, beta_1, beta_2, beta_3;
    
    #ifdef __INTEL_COMPILER
    #pragma forceinline
    #endif
    computeBeta(&beta_0, &beta_1, &beta_2, &beta_3, U_array, idx_side);
    
    /*
     * Compute the weights omega_upwind.
     */
    
    double omega_upwind_0, omega_upwind_1, omega_upwind_2;
    
    double tau_5 = fabs(beta_0 - beta_2);
    
    omega_upwind_0 = 1.0/16.0*(1.0 + pow(tau_5/(beta_0 + EPSILON), p));
    omega_upwind_1 = 5.0/8.0*(1.0 + pow(tau_5/(beta_1 + EPSILON), p));
    omega_upwind_2 = 5.0/16.0*(1.0 + pow(tau_5/(beta_2 + EPSILON), p));
    
    double omega_upwind_sum = omega_upwind_0 + omega_upwind_1 + omega_upwind_2;
    
    omega_upwind_0 = omega_upwind_0/omega_upwind_sum;
    omega_upwind_1 = omega_upwind_1/omega_upwind_sum;
    omega_upwind_2 = omega_upwind_2/omega_upwind_sum;
    
    /*
     * Compute the weights omega_central (store in omega first).
     */
    
    double omega_0, omega_1, omega_2, omega_3;
    
    double beta_avg = 1.0/8*(beta_0 + beta_2 + 6*beta_1);
    double tau_6 = fabs(beta_3 - beta_avg);
    
    omega_0 = 1.0/32.0*(C + pow(tau_6/(beta_0 + EPSILON), q));
    omega_1 = 15.0/32.0*(C + pow(tau_6/(beta_1 + EPSILON), q));
    omega_2 = 15.0/32.0*(C + pow(tau_6/(beta_2 + EPSILON), q));
    omega_3 = 1.0/32.0*(C + pow(tau_6/(beta_3 + EPSILON), q));
    
    double omega_sum = omega_0 + omega_1 + omega_2 + omega_3;
    
    omega_0 = omega_0/omega_sum;
    omega_1 = omega_1/omega_sum;
    omega_2 = omega_2/omega_sum;
    omega_3 = omega_3/omega_sum;
    
    /*
     * Compute the weights omega.
     */
    
    double R_tau = fabs(tau_6/(beta_avg + EPSILON));
    
    if (R_tau > alpha_tau)
    {
        omega_0 = sigma*omega_upwind_0 + (1.0 - sigma)*omega_0;
        omega_1 = sigma*omega_upwind_1 + (1.0 - sigma)*omega_1;
        omega_2 = sigma*omega_upwind_2 + (1.0 - sigma)*omega_2;
        omega_3 = (1.0 - sigma)*omega_3;
    }
    
    /*
     * Compute U_minus.
     */
    
    U_minus[idx_side] = 3.0/8.0*omega_0*U_array[0][idx_side] +
        (-10.0/8.0*omega_0 - 1.0/8.0*omega_1)*U_array[1][idx_side] +
        (15.0/8.0*omega_0 + 6.0/8.0*omega_1 + 3.0/8.0*omega_2)*U_array[2][idx_side] +
        (3.0/8.0*omega_1 + 6.0/8.0*omega_2 + 15.0/8.0*omega_3)*U_array[3][idx_side] +
        (-1.0/8.0*omega_2 - 10.0/8.0*omega_3)*U_array[4][idx_side] +
        3.0/8.0*omega_3*U_array[5][idx_side];
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::performWENOInterpolation_plus(
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
    const double& alpha_tau  = d_constant_alpha_tau;
    
    /*
     * Compute sigma.
     */
    
    double sigma;
    
    #ifdef __INTEL_COMPILER
    #pragma forceinline
    #endif
    computeSigma(&sigma, U_array, idx_side);
    
    /*
     * Compute beta_tilde's.
     */
    
    double beta_tilde_0, beta_tilde_1, beta_tilde_2, beta_tilde_3;
    
    #ifdef __INTEL_COMPILER
    #pragma forceinline
    #endif
    computeBetaTilde(&beta_tilde_0, &beta_tilde_1, &beta_tilde_2, &beta_tilde_3, U_array, idx_side);
    
    /*
     * Compute the weights omega_upwind_tilde.
     */
    
    double omega_upwind_tilde_0, omega_upwind_tilde_1, omega_upwind_tilde_2;
    
    double tau_5_tilde = fabs(beta_tilde_0 - beta_tilde_2);
    
    omega_upwind_tilde_0 = 1.0/16.0*(1.0 + pow(tau_5_tilde/(beta_tilde_0 + EPSILON), p));
    omega_upwind_tilde_1 = 5.0/8.0*(1.0 + pow(tau_5_tilde/(beta_tilde_1 + EPSILON), p));
    omega_upwind_tilde_2 = 5.0/16.0*(1.0 + pow(tau_5_tilde/(beta_tilde_2 + EPSILON), p));
    
    double omega_upwind_tilde_sum = omega_upwind_tilde_0 + omega_upwind_tilde_1 + omega_upwind_tilde_2;
    
    omega_upwind_tilde_0 = omega_upwind_tilde_0/omega_upwind_tilde_sum;
    omega_upwind_tilde_1 = omega_upwind_tilde_1/omega_upwind_tilde_sum;
    omega_upwind_tilde_2 = omega_upwind_tilde_2/omega_upwind_tilde_sum;
    
    /*
     * Compute the weights omega_central_tilde  (store in omega_tilde first).
     */
    
    double omega_tilde_0, omega_tilde_1, omega_tilde_2, omega_tilde_3;
    
    double beta_avg_tilde = 1.0/8*(beta_tilde_0 + beta_tilde_2 + 6*beta_tilde_1);
    double tau_6_tilde = fabs(beta_tilde_3 - beta_avg_tilde);
    
    omega_tilde_0 = 1.0/32.0*(C + pow(tau_6_tilde/(beta_tilde_0 + EPSILON), q));
    omega_tilde_1 = 15.0/32.0*(C + pow(tau_6_tilde/(beta_tilde_1 + EPSILON), q));
    omega_tilde_2 = 15.0/32.0*(C + pow(tau_6_tilde/(beta_tilde_2 + EPSILON), q));
    omega_tilde_3 = 1.0/32.0*(C + pow(tau_6_tilde/(beta_tilde_3 + EPSILON), q));
    
    double omega_tilde_sum = omega_tilde_0 + omega_tilde_1 + omega_tilde_2 + omega_tilde_3;
    
    omega_tilde_0 = omega_tilde_0/omega_tilde_sum;
    omega_tilde_1 = omega_tilde_1/omega_tilde_sum;
    omega_tilde_2 = omega_tilde_2/omega_tilde_sum;
    omega_tilde_3 = omega_tilde_3/omega_tilde_sum;
    
    /*
     * Compute the weights omega_tilde.
     */
    
    double R_tau = fabs(tau_6_tilde/(beta_avg_tilde + EPSILON));
    
    if (R_tau > alpha_tau)
    {
        omega_tilde_0 = sigma*omega_upwind_tilde_0 + (1.0 - sigma)*omega_tilde_0;
        omega_tilde_1 = sigma*omega_upwind_tilde_1 + (1.0 - sigma)*omega_tilde_1;
        omega_tilde_2 = sigma*omega_upwind_tilde_2 + (1.0 - sigma)*omega_tilde_2;
        omega_tilde_3 = (1.0 - sigma)*omega_tilde_3;
    }
    
    U_plus[idx_side] = 3.0/8.0*omega_tilde_0*U_array[5][idx_side] +
        (-10.0/8.0*omega_tilde_0 - 1.0/8.0*omega_tilde_1)*U_array[4][idx_side] +
        (15.0/8.0*omega_tilde_0 + 6.0/8.0*omega_tilde_1 + 3.0/8.0*omega_tilde_2)*U_array[3][idx_side] +
        (3.0/8.0*omega_tilde_1 + 6.0/8.0*omega_tilde_2 + 15.0/8.0*omega_tilde_3)*U_array[2][idx_side] +
        (-1.0/8.0*omega_tilde_2 - 10.0/8.0*omega_tilde_3)*U_array[1][idx_side] +
        3.0/8.0*omega_tilde_3*U_array[0][idx_side];
}
