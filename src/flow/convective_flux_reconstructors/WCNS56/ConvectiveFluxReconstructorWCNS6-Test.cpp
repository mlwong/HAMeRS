#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS6-Test.hpp"

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
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
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
         * The numbers of ghost cells are also determined.
         */
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables =
            d_flow_model->getGlobalCellDataConservativeVariables();
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > primitive_variables =
            d_flow_model->getGlobalCellDataPrimitiveVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        num_subghosts_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> num_subghosts_primitive_var;
        num_subghosts_primitive_var.reserve(d_num_eqn);
        
        std::vector<double*> Q;
        Q.reserve(d_num_eqn);
        
        std::vector<double*> V;
        V.reserve(d_num_eqn);
        
        int count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the conservative variable vector is not in the system of equations,
                // ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                Q.push_back(conservative_variables[vi]->getPointer(di));
                num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
                
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the primitive variable vector is not in the system of equations,
                // ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                V.push_back(primitive_variables[vi]->getPointer(di));
                num_subghosts_primitive_var.push_back(primitive_variables[vi]->getGhostCellWidth());
                
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
        
        boost::shared_ptr<pdat::SideData<int> > bounded_flag_minus;
        boost::shared_ptr<pdat::SideData<int> > bounded_flag_plus;
        
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
        
        bounded_flag_minus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        bounded_flag_plus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
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
         * Peform WENO interpolation.
         */
        
        t_WENO_interpolation->start();
        
        performWENOInterpolation(
            characteristic_variables_minus,
            characteristic_variables_plus,
            characteristic_variables);
        
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
         * Declare containers to store pointers for computing mid-point fluxes.
         */
        
        std::vector<double*> V_minus;
        std::vector<double*> V_plus;
        V_minus.resize(d_num_eqn);
        V_plus.resize(d_num_eqn);
        
        int* flag_minus = nullptr;
        int* flag_plus = nullptr;
        
        // Declare and initialize containers to store the references to the
        // WENO interpolated values.
        std::vector<boost::reference_wrapper<double> > V_minus_ref;
        std::vector<boost::reference_wrapper<double> > V_plus_ref;
        V_minus_ref.reserve(d_num_eqn);
        V_plus_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the mid-point flux.
        std::vector<boost::reference_wrapper<double> > F_midpoint_ref;
        F_midpoint_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the mid-point velocity.
        std::vector<boost::reference_wrapper<double> > vel_midpoint_ref;
        vel_midpoint_ref.reserve(d_dim.getValue());
        
        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */
        
        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
            bounded_flag_minus,
            primitive_variables_minus);
        
        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
            bounded_flag_plus,
            primitive_variables_plus);
        
        /*
         * Use first order interpolation if interpolated side primitive variables in x-direction
         * are out of bounds.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
        }
        
        flag_minus = bounded_flag_minus->getPointer(0);
        flag_plus = bounded_flag_plus->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            
            #ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
            #endif
            for (int i = -1; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear indices.
                const int idx_midpoint_x = i + 1;
                const int idx_cell_L = i - 1 + num_subghosts_0_primitive_var;
                const int idx_cell_R = i + num_subghosts_0_primitive_var;
                
                if (flag_minus[idx_midpoint_x] == 0 || flag_plus[idx_midpoint_x] == 0)
                {
                    V_minus[ei][idx_midpoint_x] = V[ei][idx_cell_L];
                    V_plus[ei][idx_midpoint_x] = V[ei][idx_cell_R];
                }
            }
        }
        
        /*
         * Compute mid-point flux in the x-direction.
         */
        
        t_Riemann_solver->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
        }
        
        for (int i = -1; i < interior_dim_0 + 2; i++)
        {
            // Compute the linear index of the side.
            const int idx_midpoint_x = i + 1;
            
            // Initialize container that stores the references to the mid-point flux.
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                F_midpoint_ref.push_back(boost::ref(F_midpoint_x[ei][idx_midpoint_x]));
            }
            
            // Initialize container that stores the references to the mid-point velocity.
            for (int di = 0; di < d_dim.getValue(); di++)
            {
                vel_midpoint_ref.push_back(
                    boost::ref(velocity_intercell->getPointer(0, di)[idx_midpoint_x]));
            }
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                V_minus_ref.push_back(boost::ref(V_minus[ei][idx_midpoint_x]));
                V_plus_ref.push_back(boost::ref(V_plus[ei][idx_midpoint_x]));
            }
            
            d_flow_model->
                computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                    F_midpoint_ref,
                    vel_midpoint_ref,
                    V_minus_ref,
                    V_plus_ref,
                    DIRECTION::X_DIRECTION,
                    RIEMANN_SOLVER::HLLC);
            
            F_midpoint_ref.clear();
            vel_midpoint_ref.clear();
            V_minus_ref.clear();
            V_plus_ref.clear();
        }
        
        t_Riemann_solver->stop();
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        t_reconstruct_flux->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            double* F_face_x = convective_flux->getPointer(0, ei);
            
            #ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
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
                
                #ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
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
            #ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
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
        
        std::vector<hier::IntVector> num_subghosts_primitive_var;
        num_subghosts_primitive_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_primitive_var;
        subghostcell_dims_primitive_var.reserve(d_num_eqn);
        
        std::vector<double*> Q;
        Q.reserve(d_num_eqn);
        
        std::vector<double*> V;
        V.reserve(d_num_eqn);
        
        int count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the conservative variable vector is not in the system of equations,
                // ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                Q.push_back(conservative_variables[vi]->getPointer(di));
                num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
                subghostcell_dims_conservative_var.push_back(
                    conservative_variables[vi]->getGhostBox().numberCells());
                
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the primitive variable vector is not in the system of equations,
                // ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                V.push_back(primitive_variables[vi]->getPointer(di));
                num_subghosts_primitive_var.push_back(primitive_variables[vi]->getGhostCellWidth());
                subghostcell_dims_primitive_var.push_back(
                    primitive_variables[vi]->getGhostBox().numberCells());
                
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
        
        boost::shared_ptr<pdat::SideData<int> > bounded_flag_minus;
        boost::shared_ptr<pdat::SideData<int> > bounded_flag_plus;
        
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
        
        bounded_flag_minus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        bounded_flag_plus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
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
         * Peform WENO interpolation.
         */
        
        t_WENO_interpolation->start();
        
        performWENOInterpolation(
            characteristic_variables_minus,
            characteristic_variables_plus,
            characteristic_variables);
        
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
         * Declare containers to store pointers for computing mid-point fluxes.
         */
        
        std::vector<double*> V_minus;
        std::vector<double*> V_plus;
        V_minus.resize(d_num_eqn);
        V_plus.resize(d_num_eqn);
        
        int* flag_minus = nullptr;
        int* flag_plus = nullptr;
        
        // Declare and initialize containers to store the references to the
        // WENO interpolated values.
        std::vector<boost::reference_wrapper<double> > V_minus_ref;
        std::vector<boost::reference_wrapper<double> > V_plus_ref;
        V_minus_ref.reserve(d_num_eqn);
        V_plus_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the mid-point flux.
        std::vector<boost::reference_wrapper<double> > F_midpoint_ref;
        F_midpoint_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the mid-point velocity.
        std::vector<boost::reference_wrapper<double> > vel_midpoint_ref;
        vel_midpoint_ref.reserve(d_dim.getValue());
        
        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */
        
        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
            bounded_flag_minus,
            primitive_variables_minus);
        
        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
            bounded_flag_plus,
            primitive_variables_plus);
        
        /*
         * Use first order interpolation if interpolated side primitive variables in x-direction
         * are out of bounds.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
        }
        
        flag_minus = bounded_flag_minus->getPointer(0);
        flag_plus = bounded_flag_plus->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
            const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                #ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
                #endif
                for (int i = -1; i < interior_dim_0 + 2; i++)
                {
                    // Compute the linear indices.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    const int idx_cell_L = (i - 1 + num_subghosts_0_primitive_var) +
                        (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                    
                    const int idx_cell_R = (i + num_subghosts_0_primitive_var) +
                        (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                    
                    if (flag_minus[idx_midpoint_x] == 0 || flag_plus[idx_midpoint_x] == 0)
                    {
                        V_minus[ei][idx_midpoint_x] = V[ei][idx_cell_L];
                        V_plus[ei][idx_midpoint_x] = V[ei][idx_cell_R];
                    }
                }
            }
        }
        
        /*
         * Use first order interpolation if interpolated side primitive variables in y-direction
         * are out of bounds.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(1);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(1);
        }
        
        flag_minus = bounded_flag_minus->getPointer(1);
        flag_plus = bounded_flag_plus->getPointer(1);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
            const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
            
            for (int j = -1; j < interior_dim_1 + 2; j++)
            {
                #ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
                #endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                        (j - 1 + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                    
                    const int idx_cell_T = (i + num_subghosts_0_primitive_var) +
                        (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;
                    
                    if (flag_minus[idx_midpoint_y] == 0 || flag_plus[idx_midpoint_y] == 0)
                    {
                        V_minus[ei][idx_midpoint_y] = V[ei][idx_cell_B];
                        V_plus[ei][idx_midpoint_y] = V[ei][idx_cell_T];
                    }
                }
            }
        }
        
        /*
         * Compute mid-point flux in the x-direction.
         */
        
        t_Riemann_solver->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
        }
        
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
                    F_midpoint_ref.push_back(boost::ref(F_midpoint_x[ei][idx_midpoint_x]));
                }
                
                // Initialize container that stores the references to the mid-point velocity.
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    vel_midpoint_ref.push_back(
                        boost::ref(velocity_intercell->getPointer(0, di)[idx_midpoint_x]));
                }
                
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    V_minus_ref.push_back(boost::ref(V_minus[ei][idx_midpoint_x]));
                    V_plus_ref.push_back(boost::ref(V_plus[ei][idx_midpoint_x]));
                }
                
                // Apply the Riemann solver.
                if (s > 0.65)
                {
                    d_flow_model->
                        computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                            F_midpoint_ref,
                            vel_midpoint_ref,
                            V_minus_ref,
                            V_plus_ref,
                            DIRECTION::X_DIRECTION,
                            RIEMANN_SOLVER::HLLC_HLL);
                }
                else
                {
                    d_flow_model->
                        computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                            F_midpoint_ref,
                            vel_midpoint_ref,
                            V_minus_ref,
                            V_plus_ref,
                            DIRECTION::X_DIRECTION,
                            RIEMANN_SOLVER::HLLC);
                }
                
                F_midpoint_ref.clear();
                vel_midpoint_ref.clear();
                V_minus_ref.clear();
                V_plus_ref.clear();
            }
        }
        
        t_Riemann_solver->stop();
        
        /*
         * Compute mid-point flux in the y-direction.
         */
        
        t_Riemann_solver->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(1);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(1);
        }
        
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
                    F_midpoint_ref.push_back(boost::ref(F_midpoint_y[ei][idx_midpoint_y]));
                }
                
                // Initialize container that stores the references to the mid-point velocity.
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    vel_midpoint_ref.push_back(
                        boost::ref(velocity_intercell->getPointer(1, di)[idx_midpoint_y]));
                }
                
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    V_minus_ref.push_back(boost::ref(V_minus[ei][idx_midpoint_y]));
                    V_plus_ref.push_back(boost::ref(V_plus[ei][idx_midpoint_y]));
                }
                
                // Apply the Riemann solver.
                if (s > 0.65)
                {
                    d_flow_model->
                        computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                            F_midpoint_ref,
                            vel_midpoint_ref,
                            V_minus_ref,
                            V_plus_ref,
                            DIRECTION::Y_DIRECTION,
                            RIEMANN_SOLVER::HLLC_HLL);
                }
                else
                {
                    d_flow_model->
                        computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                            F_midpoint_ref,
                            vel_midpoint_ref,
                            V_minus_ref,
                            V_plus_ref,
                            DIRECTION::Y_DIRECTION,
                            RIEMANN_SOLVER::HLLC);
                }
                
                F_midpoint_ref.clear();
                vel_midpoint_ref.clear();
                V_minus_ref.clear();
                V_plus_ref.clear();
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
                #ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
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
                #ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
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
                    #ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
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
                #ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
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
        
        std::vector<hier::IntVector> num_subghosts_primitive_var;
        num_subghosts_primitive_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
            subghostcell_dims_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_primitive_var;
            subghostcell_dims_primitive_var.reserve(d_num_eqn);
        
        std::vector<double*> Q;
        Q.reserve(d_num_eqn);
        
        std::vector<double*> V;
        V.reserve(d_num_eqn);
        
        int count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the conservative variable vector is not in the system of equations,
                // ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                Q.push_back(conservative_variables[vi]->getPointer(di));
                num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
                subghostcell_dims_conservative_var.push_back(
                    conservative_variables[vi]->getGhostBox().numberCells());
                
                count_eqn++;
            }
        }
        
        count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the primitive variable vector is not in the system of equations,
                // ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                V.push_back(primitive_variables[vi]->getPointer(di));
                num_subghosts_primitive_var.push_back(primitive_variables[vi]->getGhostCellWidth());
                subghostcell_dims_primitive_var.push_back(
                    primitive_variables[vi]->getGhostBox().numberCells());
                
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
        
        boost::shared_ptr<pdat::SideData<int> > bounded_flag_minus;
        boost::shared_ptr<pdat::SideData<int> > bounded_flag_plus;
        
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
        
        bounded_flag_minus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        bounded_flag_plus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
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
         * Peform WENO interpolation.
         */
        
        t_WENO_interpolation->start();
        
        performWENOInterpolation(
            characteristic_variables_minus,
            characteristic_variables_plus,
            characteristic_variables);
        
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
         * Declare containers to store pointers for computing mid-point fluxes.
         */
        
        std::vector<double*> V_minus;
        std::vector<double*> V_plus;
        V_minus.resize(d_num_eqn);
        V_plus.resize(d_num_eqn);
        
        int* flag_minus = nullptr;
        int* flag_plus = nullptr;
        
        // Declare and initialize containers to store the references to the
        // WENO interpolated values.
        std::vector<boost::reference_wrapper<double> > V_minus_ref;
        std::vector<boost::reference_wrapper<double> > V_plus_ref;
        V_minus_ref.reserve(d_num_eqn);
        V_plus_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the mid-point flux.
        std::vector<boost::reference_wrapper<double> > F_midpoint_ref;
        F_midpoint_ref.reserve(d_num_eqn);
        
        // Declare container to store the references to the mid-point velocity.
        std::vector<boost::reference_wrapper<double> > vel_midpoint_ref;
        vel_midpoint_ref.reserve(d_dim.getValue());
        
        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */
        
        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
            bounded_flag_minus,
            primitive_variables_minus);
        
        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
            bounded_flag_plus,
            primitive_variables_plus);
        
        /*
         * Use first order interpolation if interpolated side primitive variables in x-direction
         * are out of bounds.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
        }
        
        flag_minus = bounded_flag_minus->getPointer(0);
        flag_plus = bounded_flag_plus->getPointer(0);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
            const int num_subghosts_2_primitive_var = num_subghosts_primitive_var[ei][2];
            const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
            const int subghostcell_dim_1_primitive_var = subghostcell_dims_primitive_var[ei][1];
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
                    #endif
                    for (int i = -1; i < interior_dim_0 + 2; i++)
                    {
                        // Compute the linear indices.
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        const int idx_cell_L = (i - 1 + num_subghosts_0_primitive_var) +
                            (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                            (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                subghostcell_dim_1_primitive_var;
                        
                        const int idx_cell_R = (i + num_subghosts_0_primitive_var) +
                            (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                            (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                subghostcell_dim_1_primitive_var;
                        
                        if (flag_minus[idx_midpoint_x] == 0 || flag_plus[idx_midpoint_x] == 0)
                        {
                            V_minus[ei][idx_midpoint_x] = V[ei][idx_cell_L];
                            V_plus[ei][idx_midpoint_x] = V[ei][idx_cell_R];
                        }
                    }
                }
            }
        }
        
        /*
         * Use first order interpolation if interpolated side primitive variables in y-direction
         * are out of bounds.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(1);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(1);
        }
        
        flag_minus = bounded_flag_minus->getPointer(1);
        flag_plus = bounded_flag_plus->getPointer(1);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
            const int num_subghosts_2_primitive_var = num_subghosts_primitive_var[ei][2];
            const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
            const int subghostcell_dim_1_primitive_var = subghostcell_dims_primitive_var[ei][1];
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -1; j < interior_dim_1 + 2; j++)
                {
                    #ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
                    #endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // // Compute the linear indices.
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                            (j - 1 + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                            (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                subghostcell_dim_1_primitive_var;
                        
                        const int idx_cell_T = (i + num_subghosts_0_primitive_var) +
                            (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                            (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                subghostcell_dim_1_primitive_var;
                        
                        if (flag_minus[idx_midpoint_y] == 0 || flag_plus[idx_midpoint_y] == 0)
                        {
                            V_minus[ei][idx_midpoint_y] = V[ei][idx_cell_B];
                            V_plus[ei][idx_midpoint_y] = V[ei][idx_cell_T];
                        }
                    }
                }
            }
        }
        
        /*
         * Use first order interpolation if interpolated side primitive variables in z-direction
         * are out of bounds.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(2);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(2);
        }
        
        flag_minus = bounded_flag_minus->getPointer(2);
        flag_plus = bounded_flag_plus->getPointer(2);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
            const int num_subghosts_2_primitive_var = num_subghosts_primitive_var[ei][2];
            const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
            const int subghostcell_dim_1_primitive_var = subghostcell_dims_primitive_var[ei][1];
            
            for (int k = -1; k < interior_dim_2 + 2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
                    #endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // // Compute the linear indices.
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                            (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                            (k - 1 + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                subghostcell_dim_1_primitive_var;
                        
                        const int idx_cell_F = (i + num_subghosts_0_primitive_var) +
                            (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var +
                            (k + num_subghosts_2_primitive_var)*subghostcell_dim_0_primitive_var*
                                subghostcell_dim_1_primitive_var;
                        
                        if (flag_minus[idx_midpoint_z] == 0 || flag_plus[idx_midpoint_z] == 0)
                        {
                            V_minus[ei][idx_midpoint_z] = V[ei][idx_cell_B];
                            V_plus[ei][idx_midpoint_z] = V[ei][idx_cell_F];
                        }
                    }
                }
            }
        }
        
        /*
         * Compute mid-point flux in the x-direction.
         */
        
        t_Riemann_solver->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
        }
        
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
                        F_midpoint_ref.push_back(boost::ref(F_midpoint_x[ei][idx_midpoint_x]));
                    }
                    
                    // Initialize container that stores the references to the mid-point velocity.
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        vel_midpoint_ref.push_back(
                            boost::ref(velocity_intercell->getPointer(0, di)[idx_midpoint_x]));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        V_minus_ref.push_back(boost::ref(V_minus[ei][idx_midpoint_x]));
                        V_plus_ref.push_back(boost::ref(V_plus[ei][idx_midpoint_x]));
                    }
                    
                    // Apply the Riemann solver.
                    if (s > 0.65)
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_midpoint_ref,
                                vel_midpoint_ref,
                                V_minus_ref,
                                V_plus_ref,
                                DIRECTION::X_DIRECTION,
                                RIEMANN_SOLVER::HLLC_HLL);
                    }
                    else
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_midpoint_ref,
                                vel_midpoint_ref,
                                V_minus_ref,
                                V_plus_ref,
                                DIRECTION::X_DIRECTION,
                                RIEMANN_SOLVER::HLLC);
                    }
                    
                    F_midpoint_ref.clear();
                    vel_midpoint_ref.clear();
                    V_minus_ref.clear();
                    V_plus_ref.clear();
                }
            }
        }
        
        t_Riemann_solver->stop();
        
        /*
         * Compute mid-point flux in the y-direction.
         */
        
        t_Riemann_solver->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(1);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(1);
        }
        
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
                        F_midpoint_ref.push_back(boost::ref(F_midpoint_y[ei][idx_midpoint_y]));
                    }
                    
                    // Initialize container that stores the references to the mid-point velocity.
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        vel_midpoint_ref.push_back(
                            boost::ref(velocity_intercell->getPointer(1, di)[idx_midpoint_y]));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        V_minus_ref.push_back(boost::ref(V_minus[ei][idx_midpoint_y]));
                        V_plus_ref.push_back(boost::ref(V_plus[ei][idx_midpoint_y]));
                    }
                    
                    // Apply the Riemann solver.
                    if (s > 0.65)
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_midpoint_ref,
                                vel_midpoint_ref,
                                V_minus_ref,
                                V_plus_ref,
                                DIRECTION::Y_DIRECTION,
                                RIEMANN_SOLVER::HLLC_HLL);
                    }
                    else
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_midpoint_ref,
                                vel_midpoint_ref,
                                V_minus_ref,
                                V_plus_ref,
                                DIRECTION::Y_DIRECTION,
                                RIEMANN_SOLVER::HLLC);
                    }
                    
                    F_midpoint_ref.clear();
                    vel_midpoint_ref.clear();
                    V_minus_ref.clear();
                    V_plus_ref.clear();
                }
            }
        }
        
        t_Riemann_solver->stop();
        
        /*
         * Compute mid-point flux in the z-direction.
         */
        
        t_Riemann_solver->start();
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(2);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(2);
        }
        
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
                        F_midpoint_ref.push_back(boost::ref(F_midpoint_z[ei][idx_midpoint_z]));
                    }
                    
                    // Initialize container that stores the references to the mid-point velocity.
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        vel_midpoint_ref.push_back(
                            boost::ref(velocity_intercell->getPointer(2, di)[idx_midpoint_z]));
                    }
                    
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        V_minus_ref.push_back(boost::ref(V_minus[ei][idx_midpoint_z]));
                        V_plus_ref.push_back(boost::ref(V_plus[ei][idx_midpoint_z]));
                    }
                    
                    // Apply the Riemann solver.
                    if (s > 0.65)
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_midpoint_ref,
                                vel_midpoint_ref,
                                V_minus_ref,
                                V_plus_ref,
                                DIRECTION::Z_DIRECTION,
                                RIEMANN_SOLVER::HLLC_HLL);
                    }
                    else
                    {
                        d_flow_model->
                            computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
                                F_midpoint_ref,
                                vel_midpoint_ref,
                                V_minus_ref,
                                V_plus_ref,
                                DIRECTION::Z_DIRECTION,
                                RIEMANN_SOLVER::HLLC);
                    }
                    
                    F_midpoint_ref.clear();
                    vel_midpoint_ref.clear();
                    V_minus_ref.clear();
                    V_plus_ref.clear();
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
                    #ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
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
                    #ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
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
                    #ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
                    #endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_z = k +
                            i*(interior_dim_2 + 1) +
                            j*(interior_dim_2 + 1)*interior_dim_0;
                        
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
                        #ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
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
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorWCNS6_Test::performWENOInterpolation(
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_minus,
    std::vector<boost::shared_ptr<pdat::SideData<double> > >& variables_plus,
    const std::vector<std::vector<boost::shared_ptr<pdat::SideData<double> > > >& variables)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(variables_minus.size()) == d_num_eqn);
    TBOX_ASSERT(static_cast<int>(variables_plus.size()) == d_num_eqn);
    
    TBOX_ASSERT(static_cast<int>(variables.size()) == 6);
#endif
    
    /*
     * Get the interior dimensions.
     */
    
    const hier::IntVector interior_dims = variables_minus[0]->getBox().numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(variables_minus[ei]->getBox().numberCells() == interior_dims);
        TBOX_ASSERT(variables_plus[ei]->getBox().numberCells() == interior_dims);
        
        TBOX_ASSERT(variables_minus[ei]->getGhostCellWidth() == hier::IntVector::getOne(d_dim));
        TBOX_ASSERT(variables_plus[ei]->getGhostCellWidth() == hier::IntVector::getOne(d_dim));
    }
    
    TBOX_ASSERT(static_cast<int>(variables.size()) == 6);
    
    for (int m = 0; m < 6; m++)
    {
        TBOX_ASSERT(static_cast<int>(variables[m].size()) == d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(variables[m][ei]->getBox().numberCells() == interior_dims);
            TBOX_ASSERT(variables[m][ei]->getGhostCellWidth() == hier::IntVector::getOne(d_dim));
        }
    }
#endif
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimension.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            double* U_L = variables_minus[ei]->getPointer(0);
            double* U_R = variables_plus[ei]->getPointer(0);
            
            #ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
            #endif
            for (int i = -1; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear index of the mid-point.
                const int idx_midpoint_x = i + 1;
                
                performLocalWENOInterpolation(
                    U_L,
                    U_R,
                    U_array.data(),
                    idx_midpoint_x,
                    d_constant_C,
                    d_constant_p,
                    d_constant_q,
                    d_constant_alpha_tau);
            }
        }
        
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            double* U_L = variables_minus[ei]->getPointer(0);
            double* U_R = variables_plus[ei]->getPointer(0);
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                #ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
                #endif
                for (int i = -1; i < interior_dim_0 + 2; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    performLocalWENOInterpolation(
                        U_L,
                        U_R,
                        U_array.data(),
                        idx_midpoint_x,
                        d_constant_C,
                        d_constant_p,
                        d_constant_q,
                        d_constant_alpha_tau);
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            double* U_B = variables_minus[ei]->getPointer(1);
            double* U_T = variables_plus[ei]->getPointer(1);
            
            for (int j = -1; j < interior_dim_1 + 2; j++)
            {
                #ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
                #endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index of the mid-point.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    performLocalWENOInterpolation(
                        U_B,
                        U_T,
                        U_array.data(),
                        idx_midpoint_y,
                        d_constant_C,
                        d_constant_p,
                        d_constant_q,
                        d_constant_alpha_tau);
                }
            }
        }
        
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        /*
         * Peform WENO interpolation in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(0));
            }
            
            double* U_L = variables_minus[ei]->getPointer(0);
            double* U_R = variables_plus[ei]->getPointer(0);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
                    #endif
                    for (int i = -1; i < interior_dim_0 + 2; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        performLocalWENOInterpolation(
                            U_L,
                            U_R,
                            U_array.data(),
                            idx_midpoint_x,
                            d_constant_C,
                            d_constant_p,
                            d_constant_q,
                            d_constant_alpha_tau);
                    }
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(1));
            }
            
            double* U_B = variables_minus[ei]->getPointer(1);
            double* U_T = variables_plus[ei]->getPointer(1);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -1; j < interior_dim_1 + 2; j++)
                {
                    #ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
                    #endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        performLocalWENOInterpolation(
                            U_B,
                            U_T,
                            U_array.data(),
                            idx_midpoint_y,
                            d_constant_C,
                            d_constant_p,
                            d_constant_q,
                            d_constant_alpha_tau);;
                    }
                }
            }
        }
        
        /*
         * Peform WENO interpolation in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            std::vector<double*> U_array;
            U_array.reserve(6);
            
            for (int m = 0; m < 6; m++)
            {
                U_array.push_back(variables[m][ei]->getPointer(2));
            }
            
            double* U_B = variables_minus[ei]->getPointer(2);
            double* U_F = variables_plus[ei]->getPointer(2);
            
            for (int k = -1; k < interior_dim_2 + 2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    #ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
                    #endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the mid-point.
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        performLocalWENOInterpolation(
                            U_B,
                            U_F,
                            U_array.data(),
                            idx_midpoint_z,
                            d_constant_C,
                            d_constant_p,
                            d_constant_q,
                            d_constant_alpha_tau);
                    }
                }
            }
        }
        
    } // if (d_dim == tbox::Dimension(3))
}
