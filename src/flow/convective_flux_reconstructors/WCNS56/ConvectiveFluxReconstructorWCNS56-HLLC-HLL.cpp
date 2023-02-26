#include "flow/convective_flux_reconstructors/WCNS56/ConvectiveFluxReconstructorWCNS56-HLLC-HLL.hpp"

#define EPSILON HAMERS_EPSILON

ConvectiveFluxReconstructorWCNS56::ConvectiveFluxReconstructorWCNS56(
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
    d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*4;
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
 * Compute the convective flux and source due to splitting of convective term on a patch.
 */
void
ConvectiveFluxReconstructorWCNS56::computeConvectiveFluxAndSourceOnPatch(
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
    d_flow_model->setupBasicUtilities();
    
    HAMERS_SHARED_PTR<FlowModelRiemannSolver> riemann_solver = d_flow_model->getFlowModelRiemannSolver();
    HAMERS_SHARED_PTR<FlowModelBasicUtilities> basic_utilities = d_flow_model->getFlowModelBasicUtilities();
    
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
    HAMERS_SHARED_PTR<pdat::SideData<Real> > velocity_midpoint;
    
    if (d_has_advective_eqn_form)
    {
        velocity_midpoint.reset(new pdat::SideData<Real>(
            interior_box, d_dim.getValue(), hier::IntVector::getOne(d_dim)));
    }
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux_midpoint(
        new pdat::SideData<Real>(interior_box, d_num_eqn, hier::IntVector::getOne(d_dim)));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux_midpoint_HLLC(
        new pdat::SideData<Real>(interior_box, d_num_eqn, hier::IntVector::getOne(d_dim)));
    
    HAMERS_SHARED_PTR<pdat::SideData<Real> > convective_flux_midpoint_HLLC_HLL;
    
    if (d_dim > tbox::Dimension(1))
    {
        convective_flux_midpoint_HLLC_HLL.reset(new pdat::SideData<Real>(
            interior_box, d_num_eqn, hier::IntVector::getOne(d_dim)));
    }
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity_derivatives;
    HAMERS_SHARED_PTR<pdat::CellData<Real> > dilatation;
    HAMERS_SHARED_PTR<pdat::CellData<Real> > vorticity_magnitude;
    HAMERS_SHARED_PTR<pdat::SideData<Real> > shock_sensor;;
    
    if (d_dim > tbox::Dimension(1))
    {
        velocity_derivatives.reset(new pdat::CellData<Real>(
            interior_box, d_dim.getValue()*d_dim.getValue(), hier::IntVector::getOne(d_dim)*2));
        
        dilatation.reset(new pdat::CellData<Real>(
            interior_box, 1, hier::IntVector::getOne(d_dim)*2));
        
        vorticity_magnitude.reset(new pdat::CellData<Real>(
            interior_box, 1, hier::IntVector::getOne(d_dim)*2));
        
        shock_sensor.reset(new pdat::SideData<Real>(
            interior_box, 1, hier::IntVector::getOne(d_dim)));
    }
    
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
        
        d_flow_model->registerDerivedVariables(num_subghosts_of_data);
        
        basic_utilities->registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
            d_num_conv_ghosts,
            AVERAGING::SIMPLE);
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Get the pointers to the velocity and convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(1);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        
        hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        
        const int num_subghosts_0_velocity = num_subghosts_velocity[0];
        const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
        
        Real* u = velocity->getPointer(0);
        
        std::vector<Real*> F_node_x;
        F_node_x.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
        }
        
        std::vector<Real*> F_midpoint_x;
        F_midpoint_x.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
        }
        
        /*
         * Get the pointers to the conservative variables and primitive variables.
         * The numbers of ghost cells are also determined.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_variables =
            d_flow_model->getCellDataOfConservativeVariables();
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > primitive_variables =
            d_flow_model->getCellDataOfPrimitiveVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        num_subghosts_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> num_subghosts_primitive_var;
        num_subghosts_primitive_var.reserve(d_num_eqn);
        
        std::vector<Real*> Q;
        Q.reserve(d_num_eqn);
        
        std::vector<Real*> V;
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
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > projection_variables;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > characteristic_variables;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > characteristic_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > characteristic_variables_plus;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > primitive_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > primitive_variables_plus;
        
        HAMERS_SHARED_PTR<pdat::SideData<int> > bounded_flag_minus;
        HAMERS_SHARED_PTR<pdat::SideData<int> > bounded_flag_plus;
        
        /*
         * Initialize temporary data containers for WENO interpolation.
         */
        
        const int num_projection_var = basic_utilities->getNumberOfProjectionVariablesForPrimitiveVariables();
        projection_variables.reserve(num_projection_var);
        
        for (int vi = 0; vi < num_projection_var; vi++)
        {
            projection_variables.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
        }
        
        characteristic_variables.resize(6);
        
        for (int m = 0; m < 6; m++)
        {
            characteristic_variables[m].reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                characteristic_variables[m].push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                    interior_box, 1, hier::IntVector::getOne(d_dim)));
            }
        }
        
        characteristic_variables_minus.reserve(d_num_eqn);
        characteristic_variables_plus.reserve(d_num_eqn);
        primitive_variables_minus.reserve(d_num_eqn);
        primitive_variables_plus.reserve(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            characteristic_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            characteristic_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            primitive_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            primitive_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
        }
        
        bounded_flag_minus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        bounded_flag_plus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        /*
         * Compute the side data of the projection variables for transformation between primitive variables and
         * characteristic variables.
         */
        
        basic_utilities->computeSideDataOfProjectionVariablesForPrimitiveVariables(
            projection_variables);
        
        /*
         * Transform primitive variables to characteristic variables.
         */
        
        for (int m = 0; m < 6; m++)
        {
            basic_utilities->computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
                characteristic_variables[m],
                primitive_variables,
                projection_variables,
                m - 3);
        }
        
        /*
         * Peform WENO interpolation.
         */
        
        performWENOInterpolation(
            characteristic_variables_minus,
            characteristic_variables_plus,
            characteristic_variables);
        
        /*
         * Transform characteristic variables back to primitive variables.
         */
        
        basic_utilities->computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            primitive_variables_minus,
            characteristic_variables_minus,
            projection_variables);
        
        basic_utilities->computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            primitive_variables_plus,
            characteristic_variables_plus,
            projection_variables);
        
        /*
         * Declare containers to store pointers for computing mid-point fluxes.
         */
        
        std::vector<Real*> V_minus;
        std::vector<Real*> V_plus;
        V_minus.resize(d_num_eqn);
        V_plus.resize(d_num_eqn);
        
        int* flag_minus = nullptr;
        int* flag_plus = nullptr;
        
        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */
        
        basic_utilities->checkSideDataOfPrimitiveVariablesBounded(
            bounded_flag_minus,
            primitive_variables_minus);
        
        basic_utilities->checkSideDataOfPrimitiveVariablesBounded(
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
            
            HAMERS_PRAGMA_SIMD
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
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux_midpoint,
                velocity_midpoint,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux_midpoint,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Real* F_face_x = convective_flux->getPointer(0, ei);
            
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                const int idx_midpoint_x = i + 1;
                const int idx_midpoint_x_L = i;
                const int idx_midpoint_x_R = i + 2;
                const int idx_node_L = i - 1 + num_subghosts_0_convective_flux_x;
                const int idx_node_R = i + num_subghosts_0_convective_flux_x;
                
                F_face_x[idx_face_x] = Real(dt)*(
                    Real(1)/Real(30)*(F_midpoint_x[ei][idx_midpoint_x_R] +
                        F_midpoint_x[ei][idx_midpoint_x_L]) -
                    Real(3)/Real(10)*(F_node_x[ei][idx_node_R] +
                        F_node_x[ei][idx_node_L]) +
                    Real(23)/Real(15)*F_midpoint_x[ei][idx_midpoint_x]);
            }
        }
        
        /*
         * Compute the source.
         */
        
        if (d_has_advective_eqn_form)
        {
            Real* u_midpoint_x = velocity_midpoint->getPointer(0, 0);
            
            for (int ei = 0; ei < d_num_eqn; ei ++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
                    const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                    
                    HAMERS_PRAGMA_SIMD
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
                        
                        S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                            (Real(3)/Real(2)*(u_midpoint_x[idx_midpoint_x_R] -
                                 u_midpoint_x[idx_midpoint_x_L]) -
                             Real(3)/Real(10)*(u[idx_cell_wghost_x_R] -
                                 u[idx_cell_wghost_x_L]) +
                             Real(1)/Real(30)*(u_midpoint_x[idx_midpoint_x_RR] -
                                 u_midpoint_x[idx_midpoint_x_LL]))/Real(dx[0]));
                    }
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
         * Get the interior dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRIMITIVE_VARIABLES", d_num_conv_ghosts));
        
        d_flow_model->registerDerivedVariables(num_subghosts_of_data);
        
        basic_utilities->registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
            d_num_conv_ghosts,
            AVERAGING::SIMPLE);
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Get the pointers to the velocity and convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(2);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
        
        hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
        hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
        
        const int num_subghosts_0_velocity = num_subghosts_velocity[0];
        const int num_subghosts_1_velocity = num_subghosts_velocity[1];
        const int subghostcell_dim_0_velocity = subghostcell_dims_velocity[0];
        
        const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
        const int num_subghosts_1_convective_flux_x = num_subghosts_convective_flux_x[1];
        const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
        
        const int num_subghosts_0_convective_flux_y = num_subghosts_convective_flux_y[0];
        const int num_subghosts_1_convective_flux_y = num_subghosts_convective_flux_y[1];
        const int subghostcell_dim_0_convective_flux_y = subghostcell_dims_convective_flux_y[0];
        
        Real* u     = velocity->getPointer(0);
        Real* v     = velocity->getPointer(1);
        Real* theta = dilatation->getPointer(0);
        Real* Omega = vorticity_magnitude->getPointer(0);
        Real* s_x   = shock_sensor->getPointer(0);
        Real* s_y   = shock_sensor->getPointer(1);
        
        std::vector<Real*> F_node_x;
        std::vector<Real*> F_node_y;
        F_node_x.reserve(d_num_eqn);
        F_node_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
            F_node_y.push_back(convective_flux_node[1]->getPointer(ei));
        }
        
        std::vector<Real*> F_midpoint_x;
        std::vector<Real*> F_midpoint_y;
        std::vector<Real*> F_midpoint_HLLC_x;
        std::vector<Real*> F_midpoint_HLLC_y;
        std::vector<Real*> F_midpoint_HLLC_HLL_x;
        std::vector<Real*> F_midpoint_HLLC_HLL_y;
        F_midpoint_x.reserve(d_num_eqn);
        F_midpoint_y.reserve(d_num_eqn);
        F_midpoint_HLLC_x.reserve(d_num_eqn);
        F_midpoint_HLLC_y.reserve(d_num_eqn);
        F_midpoint_HLLC_HLL_x.reserve(d_num_eqn);
        F_midpoint_HLLC_HLL_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
            F_midpoint_y.push_back(convective_flux_midpoint->getPointer(1, ei));
        }
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_HLLC_x.push_back(convective_flux_midpoint_HLLC->getPointer(0, ei));
            F_midpoint_HLLC_y.push_back(convective_flux_midpoint_HLLC->getPointer(1, ei));
        }
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_HLLC_HLL_x.push_back(convective_flux_midpoint_HLLC_HLL->getPointer(0, ei));
            F_midpoint_HLLC_HLL_y.push_back(convective_flux_midpoint_HLLC_HLL->getPointer(1, ei));
        }
        
        /*
         * Compute the derivatives of velocity, dilatation and vorticity magnitude.
         */
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder("first order derivative in x-direction", d_dim, DIRECTION::X_DIRECTION, 1));
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
            new DerivativeFirstOrder("first order derivative in y-direction", d_dim, DIRECTION::Y_DIRECTION, 1));
        
        // Compute dudx.
        derivative_first_order_x->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[0]),
            0,
            0);
        
        // Compute dudy.
        derivative_first_order_y->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[1]),
            1,
            0);
        
        // Compute dvdx.
        derivative_first_order_x->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[0]),
            2,
            1);
        
        // Compute dvdy.
        derivative_first_order_y->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[1]),
            3,
            1);
        
        // Get the pointers to the cell data of velocity derivatives.
        Real* dudx = velocity_derivatives->getPointer(0);
        Real* dudy = velocity_derivatives->getPointer(1);
        Real* dvdx = velocity_derivatives->getPointer(2);
        Real* dvdy = velocity_derivatives->getPointer(3);
        
        // Compute the dilatation.
        for (int j = -2; j < interior_dim_1 + 2; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -2; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear index.
                const int idx = (i + 2) +
                    (j + 2)*(interior_dim_0 + 4);
                
                theta[idx] = dudx[idx] + dvdy[idx];
            }
        }
        
        // Compute the magnitude of vorticity.
        for (int j = -2; j < interior_dim_1 + 2; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -2; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear index.
                const int idx = (i + 2) +
                    (j + 2)*(interior_dim_0 + 4);
                
                Omega[idx] = std::abs(dvdx[idx] - dudy[idx]);
            }
        }
        
        /*
         * Get the pointers to the conservative variables and primitive variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_variables =
            d_flow_model->getCellDataOfConservativeVariables();
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > primitive_variables =
            d_flow_model->getCellDataOfPrimitiveVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        num_subghosts_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> num_subghosts_primitive_var;
        num_subghosts_primitive_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_primitive_var;
        subghostcell_dims_primitive_var.reserve(d_num_eqn);
        
        std::vector<Real*> Q;
        Q.reserve(d_num_eqn);
        
        std::vector<Real*> V;
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
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > projection_variables;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > characteristic_variables;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > characteristic_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > characteristic_variables_plus;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > primitive_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > primitive_variables_plus;
        
        HAMERS_SHARED_PTR<pdat::SideData<int> > bounded_flag_minus;
        HAMERS_SHARED_PTR<pdat::SideData<int> > bounded_flag_plus;
        
        /*
         * Initialize temporary data containers for WENO interpolation.
         */
        
        const int num_projection_var = basic_utilities->getNumberOfProjectionVariablesForPrimitiveVariables();
        projection_variables.reserve(num_projection_var);
        
        for (int vi = 0; vi < num_projection_var; vi++)
        {
            projection_variables.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
        }
        
        characteristic_variables.resize(6);
        
        for (int m = 0; m < 6; m++)
        {
            characteristic_variables[m].reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                characteristic_variables[m].push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                    interior_box, 1, hier::IntVector::getOne(d_dim)));
            }
        }
        
        characteristic_variables_minus.reserve(d_num_eqn);
        characteristic_variables_plus.reserve(d_num_eqn);
        primitive_variables_minus.reserve(d_num_eqn);
        primitive_variables_plus.reserve(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            characteristic_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            characteristic_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            primitive_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            primitive_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
        }
        
        bounded_flag_minus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        bounded_flag_plus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        /*
         * Compute the side data of the projection variables for transformation between primitive variables and
         * characteristic variables.
         */
        
        basic_utilities->computeSideDataOfProjectionVariablesForPrimitiveVariables(
            projection_variables);
        
        /*
         * Transform primitive variables to characteristic variables.
         */
        
        for (int m = 0; m < 6; m++)
        {
            basic_utilities->computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
                characteristic_variables[m],
                primitive_variables,
                projection_variables,
                m - 3);
        }
        
        /*
         * Peform WENO interpolation.
         */
        
        performWENOInterpolation(
            characteristic_variables_minus,
            characteristic_variables_plus,
            characteristic_variables);
        
        /*
         * Transform characteristic variables back to primitive variables.
         */
        
        basic_utilities->computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            primitive_variables_minus,
            characteristic_variables_minus,
            projection_variables);
        
        basic_utilities->computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            primitive_variables_plus,
            characteristic_variables_plus,
            projection_variables);
        
        /*
         * Declare containers to store pointers for computing mid-point fluxes.
         */
        
        std::vector<Real*> V_minus;
        std::vector<Real*> V_plus;
        V_minus.resize(d_num_eqn);
        V_plus.resize(d_num_eqn);
        
        int* flag_minus = nullptr;
        int* flag_plus = nullptr;
        
        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */
        
        basic_utilities->checkSideDataOfPrimitiveVariablesBounded(
            bounded_flag_minus,
            primitive_variables_minus);
        
        basic_utilities->checkSideDataOfPrimitiveVariablesBounded(
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
                HAMERS_PRAGMA_SIMD
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
                HAMERS_PRAGMA_SIMD
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
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                velocity_midpoint,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
            convective_flux_midpoint_HLLC_HLL,
            primitive_variables_minus,
            primitive_variables_plus,
            DIRECTION::X_DIRECTION,
            RIEMANN_SOLVER::HLLC_HLL);
        
        // Compute the Ducros-like shock sensor.
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = -1; i < interior_dim_0 + 2; i++)
            {
                // Compute the linear index of the side.
                const int idx_midpoint_x = (i + 1) +
                    (j + 1)*(interior_dim_0 + 3);
                
                const int idx_L = (i + 1) +
                    (j + 2)*(interior_dim_0 + 4);
                
                const int idx_R = (i + 2) +
                    (j + 2)*(interior_dim_0 + 4);
                
                Real theta_avg = Real(1)/Real(2)*(theta[idx_L] + theta[idx_R]);
                Real Omega_avg = Real(1)/Real(2)*(Omega[idx_L] + Omega[idx_R]);
                
                s_x[idx_midpoint_x] = -theta_avg/(std::abs(theta_avg) + Omega_avg + EPSILON);
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -1; i < interior_dim_0 + 2; i++)
                {
                    // Compute the linear index of the side.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3);
                    
                    if (s_x[idx_midpoint_x] > 0.65)
                    {
                        F_midpoint_x[ei][idx_midpoint_x] = F_midpoint_HLLC_HLL_x[ei][idx_midpoint_x];
                    }
                    else
                    {
                        F_midpoint_x[ei][idx_midpoint_x] = F_midpoint_HLLC_x[ei][idx_midpoint_x];
                    }
                }
            }
        }
        
        /*
         * Compute mid-point flux in the y-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                velocity_midpoint,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
            convective_flux_midpoint_HLLC_HLL,
            primitive_variables_minus,
            primitive_variables_plus,
            DIRECTION::Y_DIRECTION,
            RIEMANN_SOLVER::HLLC_HLL);
        
        // Compute the Ducros-like shock sensor.
        for (int j = -1; j < interior_dim_1 + 2; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute the linear index of the side.
                const int idx_midpoint_y = (i + 1) +
                    (j + 1)*(interior_dim_0 + 2);
                
                const int idx_B = (i + 2) +
                    (j + 1)*(interior_dim_0 + 4);
                
                const int idx_T = (i + 2) +
                    (j + 2)*(interior_dim_0 + 4);
                
                Real theta_avg = Real(1)/Real(2)*(theta[idx_B] + theta[idx_T]);
                Real Omega_avg = Real(1)/Real(2)*(Omega[idx_B] + Omega[idx_T]);
                
                s_y[idx_midpoint_y] = -theta_avg/(std::abs(theta_avg) + Omega_avg + EPSILON);
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int j = -1; j < interior_dim_1 + 2; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index of the side.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2);
                    
                    if (s_y[idx_midpoint_y] > 0.65)
                    {
                        F_midpoint_y[ei][idx_midpoint_y] = F_midpoint_HLLC_HLL_y[ei][idx_midpoint_y];
                    }
                    else
                    {
                        F_midpoint_y[ei][idx_midpoint_y] = F_midpoint_HLLC_y[ei][idx_midpoint_y];
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Real* F_face_x = convective_flux->getPointer(0, ei);
            
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
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
                    
                    F_face_x[idx_face_x] = Real(dt)*(
                        Real(1)/Real(30)*(F_midpoint_x[ei][idx_midpoint_x_R] +
                            F_midpoint_x[ei][idx_midpoint_x_L]) -
                        Real(3)/Real(10)*(F_node_x[ei][idx_node_R] +
                            F_node_x[ei][idx_node_L]) +
                        Real(23)/Real(15)*F_midpoint_x[ei][idx_midpoint_x]);
                }
            }
        }
        
        /*
         * Reconstruct the flux in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Real* F_face_y = convective_flux->getPointer(1, ei);
            
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
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
                    
                    F_face_y[idx_face_y] = Real(dt)*(
                        Real(1)/Real(30)*(F_midpoint_y[ei][idx_midpoint_y_T] +
                            F_midpoint_y[ei][idx_midpoint_y_B]) -
                        Real(3)/Real(10)*(F_node_y[ei][idx_node_T] +
                            F_node_y[ei][idx_node_B]) +
                        Real(23)/Real(15)*F_midpoint_y[ei][idx_midpoint_y]);
                }
            }
        }
        
        /*
         * Compute the source.
         */
        
        if (d_has_advective_eqn_form)
        {
            Real* u_midpoint_x = velocity_midpoint->getPointer(0, 0);
            Real* v_midpoint_y = velocity_midpoint->getPointer(1, 1);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
                    const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                    const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
                    const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
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
                            
                            S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                (Real(3)/Real(2)*(u_midpoint_x[idx_midpoint_x_R] -
                                     u_midpoint_x[idx_midpoint_x_L]) -
                                 Real(3)/Real(10)*(u[idx_cell_wghost_x_R] -
                                     u[idx_cell_wghost_x_L]) +
                                 Real(1)/Real(30)*(u_midpoint_x[idx_midpoint_x_RR] -
                                     u_midpoint_x[idx_midpoint_x_LL]))/Real(dx[0]) +
                                (Real(3)/Real(2)*(v_midpoint_y[idx_midpoint_y_T] -
                                     v_midpoint_y[idx_midpoint_y_B]) -
                                 Real(3)/Real(10)*(v[idx_cell_wghost_y_T] -
                                     v[idx_cell_wghost_y_B]) +
                                 Real(1)/Real(30)*(v_midpoint_y[idx_midpoint_y_TT] -
                                     v_midpoint_y[idx_midpoint_y_BB]))/Real(dx[1]));
                        }
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
         * Get the interior dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Z", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRIMITIVE_VARIABLES", d_num_conv_ghosts));
        
        d_flow_model->registerDerivedVariables(num_subghosts_of_data);
        
        basic_utilities->registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
            d_num_conv_ghosts,
            AVERAGING::SIMPLE);
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Get the pointers to the velocity and convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity = d_flow_model->getCellData("VELOCITY");
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(3);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
        convective_flux_node[2] = d_flow_model->getCellData("CONVECTIVE_FLUX_Z");
        
        hier::IntVector num_subghosts_velocity = velocity->getGhostCellWidth();
        hier::IntVector subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
        
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
        
        Real* u     = velocity->getPointer(0);
        Real* v     = velocity->getPointer(1);
        Real* w     = velocity->getPointer(2);
        Real* theta = dilatation->getPointer(0);
        Real* Omega = vorticity_magnitude->getPointer(0);
        Real* s_x   = shock_sensor->getPointer(0);
        Real* s_y   = shock_sensor->getPointer(1);
        Real* s_z   = shock_sensor->getPointer(2);
        
        std::vector<Real*> F_node_x;
        std::vector<Real*> F_node_y;
        std::vector<Real*> F_node_z;
        F_node_x.reserve(d_num_eqn);
        F_node_y.reserve(d_num_eqn);
        F_node_z.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
            F_node_y.push_back(convective_flux_node[1]->getPointer(ei));
            F_node_z.push_back(convective_flux_node[2]->getPointer(ei));
        }
        
        std::vector<Real*> F_midpoint_x;
        std::vector<Real*> F_midpoint_y;
        std::vector<Real*> F_midpoint_z;
        std::vector<Real*> F_midpoint_HLLC_x;
        std::vector<Real*> F_midpoint_HLLC_y;
        std::vector<Real*> F_midpoint_HLLC_z;
        std::vector<Real*> F_midpoint_HLLC_HLL_x;
        std::vector<Real*> F_midpoint_HLLC_HLL_y;
        std::vector<Real*> F_midpoint_HLLC_HLL_z;
        F_midpoint_x.reserve(d_num_eqn);
        F_midpoint_y.reserve(d_num_eqn);
        F_midpoint_z.reserve(d_num_eqn);
        F_midpoint_HLLC_x.reserve(d_num_eqn);
        F_midpoint_HLLC_y.reserve(d_num_eqn);
        F_midpoint_HLLC_z.reserve(d_num_eqn);
        F_midpoint_HLLC_HLL_x.reserve(d_num_eqn);
        F_midpoint_HLLC_HLL_y.reserve(d_num_eqn);
        F_midpoint_HLLC_HLL_z.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_x.push_back(convective_flux_midpoint->getPointer(0, ei));
            F_midpoint_y.push_back(convective_flux_midpoint->getPointer(1, ei));
            F_midpoint_z.push_back(convective_flux_midpoint->getPointer(2, ei));
        }
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_HLLC_x.push_back(convective_flux_midpoint_HLLC->getPointer(0, ei));
            F_midpoint_HLLC_y.push_back(convective_flux_midpoint_HLLC->getPointer(1, ei));
            F_midpoint_HLLC_z.push_back(convective_flux_midpoint_HLLC->getPointer(2, ei));
        }
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_midpoint_HLLC_HLL_x.push_back(convective_flux_midpoint_HLLC_HLL->getPointer(0, ei));
            F_midpoint_HLLC_HLL_y.push_back(convective_flux_midpoint_HLLC_HLL->getPointer(1, ei));
            F_midpoint_HLLC_HLL_z.push_back(convective_flux_midpoint_HLLC_HLL->getPointer(2, ei));
        }
        
        /*
         * Compute the derivatives of velocity, dilatation and vorticity magnitude.
         */
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
            new DerivativeFirstOrder("first order derivative in x-direction", d_dim, DIRECTION::X_DIRECTION, 1));
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
            new DerivativeFirstOrder("first order derivative in y-direction", d_dim, DIRECTION::Y_DIRECTION, 1));
        
        HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_z(
            new DerivativeFirstOrder("first order derivative in z-direction", d_dim, DIRECTION::Z_DIRECTION, 1));
        
        // Compute dudx.
        derivative_first_order_x->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[0]),
            0,
            0);
        
        // Compute dudy.
        derivative_first_order_y->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[1]),
            1,
            0);
        
        // Compute dudz.
        derivative_first_order_z->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[2]),
            2,
            0);
        
        // Compute dvdx.
        derivative_first_order_x->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[0]),
            3,
            1);
        
        // Compute dvdy.
        derivative_first_order_y->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[1]),
            4,
            1);
        
        // Compute dvdz.
        derivative_first_order_z->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[2]),
            5,
            1);
        
        // Compute dwdx.
        derivative_first_order_x->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[0]),
            6,
            2);
        
        // Compute dwdy.
        derivative_first_order_y->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[1]),
            7,
            2);
        
        // Compute dwdz.
        derivative_first_order_z->computeDerivative(
            velocity_derivatives,
            velocity,
            Real(dx[2]),
            8,
            2);
        
        // Get the pointers to the cell data of velocity derivatives.
        Real* dudx = velocity_derivatives->getPointer(0);
        Real* dudy = velocity_derivatives->getPointer(1);
        Real* dudz = velocity_derivatives->getPointer(2);
        Real* dvdx = velocity_derivatives->getPointer(3);
        Real* dvdy = velocity_derivatives->getPointer(4);
        Real* dvdz = velocity_derivatives->getPointer(5);
        Real* dwdx = velocity_derivatives->getPointer(6);
        Real* dwdy = velocity_derivatives->getPointer(7);
        Real* dwdz = velocity_derivatives->getPointer(8);
        
        // Compute the dilatation.
        for (int k = -2; k < interior_dim_2 + 2; k++)
        {
            for (int j = -2; j < interior_dim_1 + 2; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -2; i < interior_dim_0 + 2; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + 2) +
                        (j + 2)*(interior_dim_0 + 4) +
                        (k + 2)*(interior_dim_0 + 4)*
                            (interior_dim_1 + 4);
                    
                    theta[idx] = dudx[idx] + dvdy[idx] + dwdz[idx];
                }
            }
        }
        
        // Compute the magnitude of vorticity.
        for (int k = -2; k < interior_dim_2 + 2; k++)
        {
            for (int j = -2; j < interior_dim_1 + 2; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -2; i < interior_dim_0 + 2; i++)
                {
                    // Compute the linear index.
                    const int idx = (i + 2) +
                        (j + 2)*(interior_dim_0 + 4) +
                        (k + 2)*(interior_dim_0 + 4)*
                            (interior_dim_1 + 4);
                    
                    const Real omega_x = dwdy[idx] - dvdz[idx];
                    const Real omega_y = dudz[idx] - dwdx[idx];
                    const Real omega_z = dvdx[idx] - dudy[idx];
                    
                    Omega[idx] = std::sqrt(omega_x*omega_x + omega_y*omega_y + omega_z*omega_z);
                }
            }
        }
        
        /*
         * Get the pointers to the conservative variables and primitive variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > conservative_variables =
            d_flow_model->getCellDataOfConservativeVariables();
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > primitive_variables =
            d_flow_model->getCellDataOfPrimitiveVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        num_subghosts_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> num_subghosts_primitive_var;
        num_subghosts_primitive_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_primitive_var;
        subghostcell_dims_primitive_var.reserve(d_num_eqn);
        
        std::vector<Real*> Q;
        Q.reserve(d_num_eqn);
        
        std::vector<Real*> V;
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
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > projection_variables;
        
        std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > > characteristic_variables;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > characteristic_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > characteristic_variables_plus;
        
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > primitive_variables_minus;
        std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > primitive_variables_plus;
        
        HAMERS_SHARED_PTR<pdat::SideData<int> > bounded_flag_minus;
        HAMERS_SHARED_PTR<pdat::SideData<int> > bounded_flag_plus;
        
        /*
         * Initialize temporary data containers for WENO interpolation.
         */
        
        const int num_projection_var = basic_utilities->getNumberOfProjectionVariablesForPrimitiveVariables();
        projection_variables.reserve(num_projection_var);
        
        for (int vi = 0; vi < num_projection_var; vi++)
        {
            projection_variables.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
        }
        
        characteristic_variables.resize(6);
        
        for (int m = 0; m < 6; m++)
        {
            characteristic_variables[m].reserve(d_num_eqn);
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                characteristic_variables[m].push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                    interior_box, 1, hier::IntVector::getOne(d_dim)));
            }
        }
        
        characteristic_variables_minus.reserve(d_num_eqn);
        characteristic_variables_plus.reserve(d_num_eqn);
        primitive_variables_minus.reserve(d_num_eqn);
        primitive_variables_plus.reserve(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            characteristic_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            characteristic_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            primitive_variables_minus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
            
            primitive_variables_plus.push_back(HAMERS_MAKE_SHARED<pdat::SideData<Real> >(
                interior_box, 1, hier::IntVector::getOne(d_dim)));
        }
        
        bounded_flag_minus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        bounded_flag_plus.reset(
            new pdat::SideData<int>(interior_box, 1, hier::IntVector::getOne(d_dim)));
        
        /*
         * Compute the side data of the projection variables for transformation between primitive variables and
         * characteristic variables.
         */
        
        basic_utilities->computeSideDataOfProjectionVariablesForPrimitiveVariables(
            projection_variables);
        
        /*
         * Transform primitive variables to characteristic variables.
         */
        
        for (int m = 0; m < 6; m++)
        {
            basic_utilities->computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
                characteristic_variables[m],
                primitive_variables,
                projection_variables,
                m - 3);
        }
        
        /*
         * Peform WENO interpolation.
         */
        
        performWENOInterpolation(
            characteristic_variables_minus,
            characteristic_variables_plus,
            characteristic_variables);
        
        /*
         * Transform characteristic variables back to primitive variables.
         */
        
        basic_utilities->computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            primitive_variables_minus,
            characteristic_variables_minus,
            projection_variables);
        
        basic_utilities->computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            primitive_variables_plus,
            characteristic_variables_plus,
            projection_variables);
        
        /*
         * Declare containers to store pointers for computing mid-point fluxes.
         */
        
        std::vector<Real*> V_minus;
        std::vector<Real*> V_plus;
        V_minus.resize(d_num_eqn);
        V_plus.resize(d_num_eqn);
        
        int* flag_minus = nullptr;
        int* flag_plus = nullptr;
        
        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */
        
        basic_utilities->checkSideDataOfPrimitiveVariablesBounded(
            bounded_flag_minus,
            primitive_variables_minus);
        
        basic_utilities->checkSideDataOfPrimitiveVariablesBounded(
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
                    HAMERS_PRAGMA_SIMD
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
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
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
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
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
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                velocity_midpoint,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
            convective_flux_midpoint_HLLC_HLL,
            primitive_variables_minus,
            primitive_variables_plus,
            DIRECTION::X_DIRECTION,
            RIEMANN_SOLVER::HLLC_HLL);
        
        // Compute the Ducros-like shock sensor.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = -1; i < interior_dim_0 + 2; i++)
                {
                    // Compute the linear index of the side.
                    const int idx_midpoint_x = (i + 1) +
                        (j + 1)*(interior_dim_0 + 3) +
                        (k + 1)*(interior_dim_0 + 3)*
                            (interior_dim_1 + 2);
                    
                    const int idx_L = (i + 1) +
                        (j + 2)*(interior_dim_0 + 4) +
                        (k + 2)*(interior_dim_0 + 4)*
                            (interior_dim_1 + 4);
                    
                    const int idx_R = (i + 2) +
                        (j + 2)*(interior_dim_0 + 4) +
                        (k + 2)*(interior_dim_0 + 4)*
                            (interior_dim_1 + 4);
                    
                    Real theta_avg = Real(1)/Real(2)*(theta[idx_L] + theta[idx_R]);
                    Real Omega_avg = Real(1)/Real(2)*(Omega[idx_L] + Omega[idx_R]);
                    
                    s_x[idx_midpoint_x] = -theta_avg/(std::abs(theta_avg) + Omega_avg + EPSILON);
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = -1; i < interior_dim_0 + 2; i++)
                    {
                        // Compute the linear index of the side.
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*
                                (interior_dim_1 + 2);
                        
                        if (s_x[idx_midpoint_x] > 0.65)
                        {
                            F_midpoint_x[ei][idx_midpoint_x] = F_midpoint_HLLC_HLL_x[ei][idx_midpoint_x];
                        }
                        else
                        {
                            F_midpoint_x[ei][idx_midpoint_x] = F_midpoint_HLLC_x[ei][idx_midpoint_x];
                        }
                    }
                }
            }
        }
        
        /*
         * Compute mid-point flux in the y-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                velocity_midpoint,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
            convective_flux_midpoint_HLLC_HLL,
            primitive_variables_minus,
            primitive_variables_plus,
            DIRECTION::Y_DIRECTION,
            RIEMANN_SOLVER::HLLC_HLL);
        
        // Compute the Ducros-like shock sensor.
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = -1; j < interior_dim_1 + 2; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index of the side.
                    const int idx_midpoint_y = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2) +
                        (k + 1)*(interior_dim_0 + 2)*
                            (interior_dim_1 + 3);
                    
                    const int idx_B = (i + 2) +
                        (j + 1)*(interior_dim_0 + 4) +
                        (k + 2)*(interior_dim_0 + 4)*
                            (interior_dim_1 + 4);
                    
                    const int idx_T = (i + 2) +
                        (j + 2)*(interior_dim_0 + 4) +
                        (k + 2)*(interior_dim_0 + 4)*
                            (interior_dim_1 + 4);
                    
                    Real theta_avg = Real(1)/Real(2)*(theta[idx_B] + theta[idx_T]);
                    Real Omega_avg = Real(1)/Real(2)*(Omega[idx_B] + Omega[idx_T]);
                    
                    s_y[idx_midpoint_y] = -theta_avg/(std::abs(theta_avg) + Omega_avg + EPSILON);
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = -1; j < interior_dim_1 + 2; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the side.
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 3);
                        
                        if (s_y[idx_midpoint_y] > 0.65)
                        {
                            F_midpoint_y[ei][idx_midpoint_y] = F_midpoint_HLLC_HLL_y[ei][idx_midpoint_y];
                        }
                        else
                        {
                            F_midpoint_y[ei][idx_midpoint_y] = F_midpoint_HLLC_y[ei][idx_midpoint_y];
                        }
                    }
                }
            }
        }
        
        /*
         * Compute mid-point flux in the z-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                velocity_midpoint,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Z_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux_midpoint_HLLC,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Z_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
            convective_flux_midpoint_HLLC_HLL,
            primitive_variables_minus,
            primitive_variables_plus,
            DIRECTION::Z_DIRECTION,
            RIEMANN_SOLVER::HLLC_HLL);
        
        // Compute the Ducros-like shock sensor.
        for (int k = -1; k < interior_dim_2 + 2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear index of the side.
                    const int idx_midpoint_z = (i + 1) +
                        (j + 1)*(interior_dim_0 + 2) +
                        (k + 1)*(interior_dim_0 + 2)*
                            (interior_dim_1 + 2);
                    
                    const int idx_B = (i + 2) +
                        (j + 2)*(interior_dim_0 + 4) +
                        (k + 1)*(interior_dim_0 + 4)*
                            (interior_dim_1 + 4);
                    
                    const int idx_F = (i + 2) +
                        (j + 2)*(interior_dim_0 + 4) +
                        (k + 2)*(interior_dim_0 + 4)*
                            (interior_dim_1 + 4);
                    
                    Real theta_avg = Real(1)/Real(2)*(theta[idx_B] + theta[idx_F]);
                    Real Omega_avg = Real(1)/Real(2)*(Omega[idx_B] + Omega[idx_F]);
                    
                    s_z[idx_midpoint_z] = -theta_avg/(std::abs(theta_avg) + Omega_avg + EPSILON);
                }
            }
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int k = -1; k < interior_dim_2 + 2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear index of the side.
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*
                                (interior_dim_1 + 2);
                        
                        if (s_z[idx_midpoint_z] > 0.65)
                        {
                            F_midpoint_z[ei][idx_midpoint_z] = F_midpoint_HLLC_HLL_z[ei][idx_midpoint_z];
                        }
                        else
                        {
                            F_midpoint_z[ei][idx_midpoint_z] = F_midpoint_HLLC_z[ei][idx_midpoint_z];
                        }
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Real* F_face_x = convective_flux->getPointer(0, ei);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dim_0 + 1) +
                            k*(interior_dim_0 + 1)*interior_dim_1;
                        
                        const int idx_midpoint_x = (i + 1) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*(interior_dim_1 + 2);
                        
                        const int idx_midpoint_x_L = i +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*(interior_dim_1 + 2);
                        
                        const int idx_midpoint_x_R = (i + 2) +
                            (j + 1)*(interior_dim_0 + 3) +
                            (k + 1)*(interior_dim_0 + 3)*(interior_dim_1 + 2);
                        
                        const int idx_node_L = (i - 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                            (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                subghostcell_dim_1_convective_flux_x;
                        
                        const int idx_node_R = (i + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                            (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                subghostcell_dim_1_convective_flux_x;
                        
                        F_face_x[idx_face_x] = Real(dt)*(
                            Real(1)/Real(30)*(F_midpoint_x[ei][idx_midpoint_x_R] +
                                F_midpoint_x[ei][idx_midpoint_x_L]) -
                            Real(3)/Real(10)*(F_node_x[ei][idx_node_R] +
                                F_node_x[ei][idx_node_L]) +
                            Real(23)/Real(15)*F_midpoint_x[ei][idx_midpoint_x]);
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in the y-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Real* F_face_y = convective_flux->getPointer(1, ei);
            
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_midpoint_y = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*(interior_dim_1 + 3);
                        
                        const int idx_midpoint_y_B = (i + 1) +
                            j*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*(interior_dim_1 + 3);
                        
                        const int idx_midpoint_y_T = (i + 1) +
                            (j + 2)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*(interior_dim_1 + 3);
                        
                        const int idx_node_B = (i + num_subghosts_0_convective_flux_y) +
                            (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                            (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                subghostcell_dim_1_convective_flux_y;
                        
                        const int idx_node_T = (i + num_subghosts_0_convective_flux_y) +
                            (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                            (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                subghostcell_dim_1_convective_flux_y;
                        
                        F_face_y[idx_face_y] = Real(dt)*(
                            Real(1)/Real(30)*(F_midpoint_y[ei][idx_midpoint_y_T] +
                                F_midpoint_y[ei][idx_midpoint_y_B]) -
                            Real(3)/Real(10)*(F_node_y[ei][idx_node_T] +
                                F_node_y[ei][idx_node_B]) +
                            Real(23)/Real(15)*F_midpoint_y[ei][idx_midpoint_y]);
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in the z-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Real* F_face_z = convective_flux->getPointer(2, ei);
            
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_z = i +
                            j*interior_dim_0 +
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_midpoint_z = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 1)*(interior_dim_0 + 2)*(interior_dim_1 + 2);
                        
                        const int idx_midpoint_z_B = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            k*(interior_dim_0 + 2)*(interior_dim_1 + 2);
                        
                        const int idx_midpoint_z_F = (i + 1) +
                            (j + 1)*(interior_dim_0 + 2) +
                            (k + 2)*(interior_dim_0 + 2)*(interior_dim_1 + 2);
                        
                        const int idx_node_B = (i + num_subghosts_0_convective_flux_z) +
                            (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                            (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                subghostcell_dim_1_convective_flux_z;
                        
                        const int idx_node_F = (i + num_subghosts_0_convective_flux_z) +
                            (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                            (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                subghostcell_dim_1_convective_flux_z;
                        
                        F_face_z[idx_face_z] = Real(dt)*(
                            Real(1)/Real(30)*(F_midpoint_z[ei][idx_midpoint_z_F] +
                                F_midpoint_z[ei][idx_midpoint_z_B]) -
                            Real(3)/Real(10)*(F_node_z[ei][idx_node_F] +
                                F_node_z[ei][idx_node_B]) +
                            Real(23)/Real(15)*F_midpoint_z[ei][idx_midpoint_z]);
                    }
                }
            }
        }
        
        /*
         * Compute the source.
         */
        
        if (d_has_advective_eqn_form)
        {
            Real* u_midpoint_x = velocity_midpoint->getPointer(0, 0);
            Real* v_midpoint_y = velocity_midpoint->getPointer(1, 1);
            Real* w_midpoint_z = velocity_midpoint->getPointer(2, 2);
            
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
                    const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
                    const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
                    const int num_subghosts_2_conservative_var = num_subghosts_conservative_var[ei][2];
                    const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
                    const int subghostcell_dim_1_conservative_var = subghostcell_dims_conservative_var[ei][1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
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
                                
                                S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                    (Real(3)/Real(2)*(u_midpoint_x[idx_midpoint_x_R] -
                                         u_midpoint_x[idx_midpoint_x_L]) -
                                     Real(3)/Real(10)*(u[idx_cell_wghost_x_R] -
                                         u[idx_cell_wghost_x_L]) +
                                     Real(1)/Real(30)*(u_midpoint_x[idx_midpoint_x_RR] -
                                         u_midpoint_x[idx_midpoint_x_LL]))/Real(dx[0]) +
                                    (Real(3)/Real(2)*(v_midpoint_y[idx_midpoint_y_T] -
                                         v_midpoint_y[idx_midpoint_y_B]) -
                                     Real(3)/Real(10)*(v[idx_cell_wghost_y_T] -
                                         v[idx_cell_wghost_y_B]) +
                                     Real(1)/Real(30)*(v_midpoint_y[idx_midpoint_y_TT] -
                                         v_midpoint_y[idx_midpoint_y_BB]))/Real(dx[1]) +
                                    (Real(3)/Real(2)*(w_midpoint_z[idx_midpoint_z_F] -
                                         w_midpoint_z[idx_midpoint_z_B]) -
                                     Real(3)/Real(10)*(w[idx_cell_wghost_z_F] -
                                         w[idx_cell_wghost_z_B]) +
                                     Real(1)/Real(30)*(w_midpoint_z[idx_midpoint_z_FF] -
                                         w_midpoint_z[idx_midpoint_z_BB]))/Real(dx[2]));
                            }
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
