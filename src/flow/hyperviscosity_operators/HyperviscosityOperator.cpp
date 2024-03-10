#include "flow/hyperviscosity_operators/HyperviscosityOperator.hpp"


HyperviscosityOperator::HyperviscosityOperator(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& hyperviscosity_operator_db):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_hyperviscosity_op_ghosts(hier::IntVector::getZero(d_dim)),
        d_num_eqn(num_eqn),
        d_flow_model(flow_model),
        d_hyperviscosity_operator_db(
            hyperviscosity_operator_db)
{
    d_lap_order = d_hyperviscosity_operator_db->
        getIntegerWithDefault("lap_order", 6);
    d_lap_order = d_hyperviscosity_operator_db->
        getIntegerWithDefault("d_lap_order", d_lap_order);
    
    d_accuracy_order = d_hyperviscosity_operator_db->
        getIntegerWithDefault("accuracy_order", 6);
    d_accuracy_order = d_hyperviscosity_operator_db->
        getIntegerWithDefault("d_accuracy_order", d_accuracy_order);
    
    d_use_flux_form = d_hyperviscosity_operator_db->
        getBoolWithDefault("use_flux_form", false);
    d_use_flux_form = d_hyperviscosity_operator_db->
        getBoolWithDefault("d_use_flux_form", d_use_flux_form);
    
    d_coeff = d_hyperviscosity_operator_db->
        getRealWithDefault("coeff", Real(1));
    d_coeff = d_hyperviscosity_operator_db->
        getRealWithDefault("d_coeff", d_coeff);
    
    if (d_accuracy_order != 2 && d_accuracy_order != 4  && d_accuracy_order != 6)
    {
        TBOX_ERROR("HyperviscosityOperator::HyperviscosityOperator:"
            " Only 2nd, 4th, 6th order accurate schemes are implemented!");
    }
    
    if (d_lap_order == 2)
    {
        if (d_accuracy_order == 2)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim);
        }
        else if (d_accuracy_order == 4)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*2;
        }
        else if (d_accuracy_order == 6)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*3;
        }
    }
    else if (d_lap_order == 4)
    {
        if (d_accuracy_order == 2)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*2;
        }
        else if (d_accuracy_order == 4)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*3;
        }
        else if (d_accuracy_order == 6)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*4;
        }
    }
    else if (d_lap_order == 6)
    {
        if (d_accuracy_order == 2)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*3;
        }
        else if (d_accuracy_order == 4)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*4;
        }
        else if (d_accuracy_order == 6)
        {
            d_num_hyperviscosity_op_ghosts = hier::IntVector::getOne(d_dim)*5;
        }
    }
    else
    {
        TBOX_ERROR("HyperviscosityOperator::HyperviscosityOperator:"
            " Only 2nd, 4th, 6th order Laplacian are implemented!");
    }
    
    d_coeffs_node.resize(6, Real(0));
    d_coeffs_midpoint.resize(5, Real(0));
    
    Real a_n = Real(0);
    Real b_n = Real(0);
    Real c_n = Real(0);
    Real d_n = Real(0);
    Real e_n = Real(0);
    Real f_n = Real(0);
    
    if (d_lap_order == 2)
    {
        if (d_accuracy_order == 2)
        {
            a_n = -Real(2);
            b_n =  Real(1);
        }
        else if (d_accuracy_order == 4)
        {
            a_n = -Real(5)/Real(2);
            b_n =  Real(4)/Real(3);
            c_n = -Real(1)/Real(12);
        }
        else if (d_accuracy_order == 6)
        {
            a_n = -Real(49)/Real(18);
            b_n =  Real(3)/Real(2);
            c_n = -Real(3)/Real(20);
            d_n =  Real(1)/Real(90);
        }
    }
    else if (d_lap_order == 4)
    {
        if (d_accuracy_order == 2)
        {
            a_n =  Real(6);
            b_n = -Real(4);
            c_n =  Real(1);
        }
        else if (d_accuracy_order == 4)
        {
            a_n =  Real(28)/Real(3);
            b_n = -Real(13)/Real(2);
            c_n =  Real(2);
            d_n = -Real(1)/Real(6);
        }
        else if (d_accuracy_order == 6)
        {
            a_n =  Real(91)/Real(8);
            b_n = -Real(122)/Real(15);
            c_n =  Real(169)/Real(60);
            d_n = -Real(2)/Real(5);
            e_n =  Real(7)/Real(240);
        }
    }
    else if (d_lap_order == 6)
    {
        if (d_accuracy_order == 2)
        {
            a_n = -Real(20);
            b_n =  Real(15);
            c_n = -Real(6);
            d_n =  Real(1);
        }
        else if (d_accuracy_order == 4)
        {
            a_n = -Real(75)/Real(2);
            b_n =  Real(29);
            c_n = -Real(13);
            d_n =  Real(3);
            e_n = -Real(1)/Real(4);
        }
        else if (d_accuracy_order == 6)
        {
            a_n = -Real(1023)/Real(20);
            b_n =  Real(323)/Real(8);
            c_n = -Real(39)/Real(2);
            d_n =  Real(87)/Real(16);
            e_n = -Real(19)/Real(24);
            f_n =  Real(13)/Real(240);
        }
    }
    
    d_coeffs_node[0] = a_n;
    d_coeffs_node[1] = b_n;
    d_coeffs_node[2] = c_n;
    d_coeffs_node[3] = d_n;
    d_coeffs_node[4] = e_n;
    d_coeffs_node[5] = f_n;
    
    const Real a_m = b_n + c_n + d_n + e_n + f_n;
    const Real b_m = c_n + d_n + e_n + f_n;
    const Real c_m = d_n + e_n + f_n;
    const Real d_m = e_n + f_n;
    const Real e_m = f_n;
    
    d_coeffs_midpoint[0] = a_m;
    d_coeffs_midpoint[1] = b_m;
    d_coeffs_midpoint[2] = c_m;
    d_coeffs_midpoint[3] = d_m;
    d_coeffs_midpoint[4] = e_m;
    
    d_prefactor = d_coeff*std::pow(-1, d_lap_order/2 + 1);
}


/*
 * Print all characteristics of the non-conservative diffusive flux divergence operator class.
 */
void
HyperviscosityOperator::printClassData(
    std::ostream& os) const
{
    os << "\nPrint HyperviscosityOperator object..."
       << std::endl;
    
    os << std::endl;
    
    os << "HyperviscosityOperator: this = "
       << (HyperviscosityOperator *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_lap_order = "
       << d_lap_order
       << std::endl;
    os << "d_accuracy_order = "
       << d_accuracy_order
       << std::endl;
    os << "d_use_flux_form = "
       << d_use_flux_form
       << std::endl;
    os << "d_coeff = "
       << d_coeff
       << std::endl;
}


/*
 * Put the characteristics of the hyperviscosity operator into the restart database.
 */
void
HyperviscosityOperator::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putInteger("d_lap_order", d_lap_order);
    restart_db->putInteger("d_accuracy_order", d_accuracy_order);
    restart_db->putReal("d_coeff", d_coeff);
}


/*
 * Perform the hyperviscosity operator on a patch.
 */
void
HyperviscosityOperator::performHyperviscosityOperatorOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<hier::CoarseFineBoundary> coarse_fine_bdry,
    const HAMERS_SHARED_PTR<pdat::SideVariable<Real> >& variable_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_source,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(variable_convective_flux);
    TBOX_ASSERT(variable_source);
#endif
    
    NULL_USE(coarse_fine_bdry);
    NULL_USE(time);
    NULL_USE(dt);
    NULL_USE(RK_step_number);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();

    // Set the domain.
    const hier::Box domain(interior_box);
    
    // Get the cell data of source term.
    HAMERS_SHARED_PTR<pdat::CellData<Real> > source(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<Real>, hier::PatchData>(
            patch.getPatchData(variable_source, data_context)));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(source);
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    if (d_use_flux_form)
    {
        TBOX_ERROR("HyperviscosityOperator::performHyperviscosityOperatorOnPatch:"
            " Flux form is not implemented yet!");
    }
    else
    {
        performHyperviscosityOperatorOnPatchSourceForm(
            patch,
            source,
            data_context,
            dt,
            domain);
    }
}


/*
 * Perform the hyperviscosity operator on a patch using source form.
 */
void
HyperviscosityOperator::performHyperviscosityOperatorOnPatchSourceForm(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::CellData<Real> > source,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double dt,
    const hier::Box& domain) const
{
    // Get the coefficients.
    const Real a_n = d_coeffs_node[0];
    const Real b_n = d_coeffs_node[1];
    const Real c_n = d_coeffs_node[2];
    const Real d_n = d_coeffs_node[3];
    const Real e_n = d_coeffs_node[4];
    const Real f_n = d_coeffs_node[5];
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
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
    
    /*
     * Get the local lower index and number of cells in each direction of the domain.
     */
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    domain_lo = domain.lower() - interior_box.lower();
    domain_dims = domain.numberCells();
    
    /*
     * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
     */
    
    d_flow_model->registerPatchWithDataContext(patch, data_context);
    
    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
    
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("SOUND_SPEED", hier::IntVector::getZero(d_dim)));
    
    d_flow_model->registerDerivedVariables(num_subghosts_of_data);
    
    d_flow_model->allocateMemoryForDerivedCellData();
    
    d_flow_model->computeDerivedCellData();
    
    HAMERS_SHARED_PTR<pdat::CellData<Real> > sound_speed = d_flow_model->getCellData("SOUND_SPEED");
    
    const hier::IntVector num_subghosts_sound_speed = sound_speed->getGhostCellWidth();
    const hier::IntVector subghostcell_dims_sound_speed = sound_speed->getGhostBox().numberCells();
    
    Real* c = sound_speed->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the local lower indices and the number of cells in each dimension.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_dim_0 = domain_dims[0];
        
        /*
         * Get the grid spacing.
         */
        const Real inv_dx_0 = Real(1)/Real(dx[0]);
        
        /*
         * Get the number of ghost cell.
         */
        const int num_subghosts_0_sound_speed = num_subghosts_sound_speed[0];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Real* S = source->getPointer(ei);
            
            const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
            
            if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*1)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                    const int idx_cell_nghost = i;
                    
                    const int idx_cell_sound_speed = i + num_subghosts_0_sound_speed;
                    const int idx_cell_wghost_x_L  = i - 1 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_R  = i + 1 + num_subghosts_0_conservative_var;
                    
                    S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                        inv_dx_0*(
                         a_n*Q[ei][idx_cell_wghost] +
                         b_n*(Q[ei][idx_cell_wghost_x_L] + Q[ei][idx_cell_wghost_x_R])
                        )
                    );
                }
            }
            else if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*2)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                    const int idx_cell_nghost = i;
                    
                    const int idx_cell_sound_speed = i + num_subghosts_0_sound_speed;
                    const int idx_cell_wghost_x_LL = i - 2 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_L  = i - 1 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_R  = i + 1 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_RR = i + 2 + num_subghosts_0_conservative_var;
                    
                    S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                        inv_dx_0*(
                         a_n*Q[ei][idx_cell_wghost] +
                         b_n*(Q[ei][idx_cell_wghost_x_L]  + Q[ei][idx_cell_wghost_x_R]) +
                         c_n*(Q[ei][idx_cell_wghost_x_LL] + Q[ei][idx_cell_wghost_x_RR])
                        )
                    );
                }
            }
            else if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*3)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                    const int idx_cell_nghost = i;
                    
                    const int idx_cell_sound_speed  = i + num_subghosts_0_sound_speed;
                    const int idx_cell_wghost_x_LLL = i - 3 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_LL  = i - 2 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_L   = i - 1 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_R   = i + 1 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_RR  = i + 2 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_RRR = i + 3 + num_subghosts_0_conservative_var;
                    
                    S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                        inv_dx_0*(
                         a_n*Q[ei][idx_cell_wghost] +
                         b_n*(Q[ei][idx_cell_wghost_x_L]   + Q[ei][idx_cell_wghost_x_R]) +
                         c_n*(Q[ei][idx_cell_wghost_x_LL]  + Q[ei][idx_cell_wghost_x_RR]) +
                         d_n*(Q[ei][idx_cell_wghost_x_LLL] + Q[ei][idx_cell_wghost_x_RRR])
                        )
                    );
                }
            }
            else if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*4)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                    const int idx_cell_nghost = i;
                    
                    const int idx_cell_sound_speed   = i + num_subghosts_0_sound_speed;
                    const int idx_cell_wghost_x_LLLL = i - 4 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_LLL  = i - 3 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_LL   = i - 2 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_L    = i - 1 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_R    = i + 1 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_RR   = i + 2 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_RRR  = i + 3 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_RRRR = i + 4 + num_subghosts_0_conservative_var;
                    
                    S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                        inv_dx_0*(
                         a_n*Q[ei][idx_cell_wghost] +
                         b_n*(Q[ei][idx_cell_wghost_x_L]    + Q[ei][idx_cell_wghost_x_R]) +
                         c_n*(Q[ei][idx_cell_wghost_x_LL]   + Q[ei][idx_cell_wghost_x_RR]) +
                         d_n*(Q[ei][idx_cell_wghost_x_LLL]  + Q[ei][idx_cell_wghost_x_RRR]) +
                         e_n*(Q[ei][idx_cell_wghost_x_LLLL] + Q[ei][idx_cell_wghost_x_RRRR])
                        )
                    );
                }
            }
            else if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*5)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_cell_wghost = i + num_subghosts_0_conservative_var;
                    const int idx_cell_nghost = i;
                    
                    const int idx_cell_sound_speed    = i + num_subghosts_0_sound_speed;
                    const int idx_cell_wghost_x_LLLLL = i - 5 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_LLLL  = i - 4 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_LLL   = i - 3 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_LL    = i - 2 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_L     = i - 1 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_R     = i + 1 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_RR    = i + 2 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_RRR   = i + 3 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_RRRR  = i + 4 + num_subghosts_0_conservative_var;
                    const int idx_cell_wghost_x_RRRRR = i + 5 + num_subghosts_0_conservative_var;
                    
                    S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                        inv_dx_0*(
                         a_n*Q[ei][idx_cell_wghost] +
                         b_n*(Q[ei][idx_cell_wghost_x_L]     + Q[ei][idx_cell_wghost_x_R]) +
                         c_n*(Q[ei][idx_cell_wghost_x_LL]    + Q[ei][idx_cell_wghost_x_RR]) +
                         d_n*(Q[ei][idx_cell_wghost_x_LLL]   + Q[ei][idx_cell_wghost_x_RRR]) +
                         e_n*(Q[ei][idx_cell_wghost_x_LLLL]  + Q[ei][idx_cell_wghost_x_RRRR]) +
                         f_n*(Q[ei][idx_cell_wghost_x_LLLLL] + Q[ei][idx_cell_wghost_x_RRRRR])
                        )
                    );
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the local lower indices and the number of cells in each dimension.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Get the grid spacings.
         */
        const Real inv_dx_0 = Real(1)/Real(dx[0]);
        const Real inv_dx_1 = Real(1)/Real(dx[1]);
        
        /*
         * Get the number of ghost cells.
         */
        const int num_subghosts_0_sound_speed = num_subghosts_sound_speed[0];
        const int num_subghosts_1_sound_speed = num_subghosts_sound_speed[1];
        const int subghostcell_dim_0_sound_speed = subghostcell_dims_sound_speed[0];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Real* S = source->getPointer(ei);
            
            const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
            const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
            const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
            
            if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*1)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_nghost = i + j*interior_dim_0;
                        
                        const int idx_cell_sound_speed = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_B = (i + num_subghosts_0_conservative_var) +
                            (j - 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_T = (i + num_subghosts_0_conservative_var) +
                            (j + 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                            inv_dx_0*(
                             a_n*Q[ei][idx_cell_wghost] +
                             b_n*(Q[ei][idx_cell_wghost_x_L] + Q[ei][idx_cell_wghost_x_R])
                            ) +
                            inv_dx_1*(
                             a_n*Q[ei][idx_cell_wghost] +
                             b_n*(Q[ei][idx_cell_wghost_y_B] + Q[ei][idx_cell_wghost_y_T])
                            )
                        );
                    }
                }
            }
            else if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*2)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_nghost = i + j*interior_dim_0;
                        
                        const int idx_cell_sound_speed = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_cell_wghost_x_LL = (i - 2 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_RR = (i + 2 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_BB = (i + num_subghosts_0_conservative_var) +
                            (j - 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_B = (i + num_subghosts_0_conservative_var) +
                            (j - 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_T = (i + num_subghosts_0_conservative_var) +
                            (j + 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_TT = (i + num_subghosts_0_conservative_var) +
                            (j + 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                            inv_dx_0*(
                             a_n*Q[ei][idx_cell_wghost] +
                             b_n*(Q[ei][idx_cell_wghost_x_L]  + Q[ei][idx_cell_wghost_x_R]) +
                             c_n*(Q[ei][idx_cell_wghost_x_LL] + Q[ei][idx_cell_wghost_x_RR])
                            ) +
                            inv_dx_1*(
                             a_n*Q[ei][idx_cell_wghost] +
                             b_n*(Q[ei][idx_cell_wghost_y_B]  + Q[ei][idx_cell_wghost_y_T]) +
                             c_n*(Q[ei][idx_cell_wghost_y_BB] + Q[ei][idx_cell_wghost_y_TT])
                            )
                        );
                    }
                }
            }
            else if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*3)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_nghost = i + j*interior_dim_0;
                        
                        const int idx_cell_sound_speed = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_cell_wghost_x_LLL = (i - 3 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_LL = (i - 2 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_RR = (i + 2 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_RRR = (i + 3 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_BBB = (i + num_subghosts_0_conservative_var) +
                            (j - 3 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_BB = (i + num_subghosts_0_conservative_var) +
                            (j - 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_B = (i + num_subghosts_0_conservative_var) +
                            (j - 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_T = (i + num_subghosts_0_conservative_var) +
                            (j + 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_TT = (i + num_subghosts_0_conservative_var) +
                            (j + 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_TTT = (i + num_subghosts_0_conservative_var) +
                            (j + 3 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                            inv_dx_0*(
                             a_n*Q[ei][idx_cell_wghost] +
                             b_n*(Q[ei][idx_cell_wghost_x_L]   + Q[ei][idx_cell_wghost_x_R]) +
                             c_n*(Q[ei][idx_cell_wghost_x_LL]  + Q[ei][idx_cell_wghost_x_RR]) +
                             d_n*(Q[ei][idx_cell_wghost_x_LLL] + Q[ei][idx_cell_wghost_x_RRR])
                            ) +
                            inv_dx_1*(
                             a_n*Q[ei][idx_cell_wghost] +
                             b_n*(Q[ei][idx_cell_wghost_y_B]   + Q[ei][idx_cell_wghost_y_T]) +
                             c_n*(Q[ei][idx_cell_wghost_y_BB]  + Q[ei][idx_cell_wghost_y_TT]) +
                             d_n*(Q[ei][idx_cell_wghost_y_BBB] + Q[ei][idx_cell_wghost_y_TTT])
                            )
                        );
                    }
                }
            }
            else if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*4)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_nghost = i + j*interior_dim_0;
                        
                        const int idx_cell_sound_speed = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_cell_wghost_x_LLLL = (i - 4 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_LLL = (i - 3 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_LL = (i - 2 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_RR = (i + 2 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_RRR = (i + 3 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_RRRR = (i + 4 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_BBBB = (i + num_subghosts_0_conservative_var) +
                            (j - 4 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_BBB = (i + num_subghosts_0_conservative_var) +
                            (j - 3 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_BB = (i + num_subghosts_0_conservative_var) +
                            (j - 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_B = (i + num_subghosts_0_conservative_var) +
                            (j - 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_T = (i + num_subghosts_0_conservative_var) +
                            (j + 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_TT = (i + num_subghosts_0_conservative_var) +
                            (j + 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_TTT = (i + num_subghosts_0_conservative_var) +
                            (j + 3 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_TTTT = (i + num_subghosts_0_conservative_var) +
                            (j + 4 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                            inv_dx_0*(
                             a_n*Q[ei][idx_cell_wghost] +
                             b_n*(Q[ei][idx_cell_wghost_x_L]    + Q[ei][idx_cell_wghost_x_R]) +
                             c_n*(Q[ei][idx_cell_wghost_x_LL]   + Q[ei][idx_cell_wghost_x_RR]) +
                             d_n*(Q[ei][idx_cell_wghost_x_LLL]  + Q[ei][idx_cell_wghost_x_RRR]) +
                             e_n*(Q[ei][idx_cell_wghost_x_LLLL] + Q[ei][idx_cell_wghost_x_RRRR])
                            ) +
                            inv_dx_1*(
                             a_n*Q[ei][idx_cell_wghost] +
                             b_n*(Q[ei][idx_cell_wghost_y_B]    + Q[ei][idx_cell_wghost_y_T]) +
                             c_n*(Q[ei][idx_cell_wghost_y_BB]   + Q[ei][idx_cell_wghost_y_TT]) +
                             d_n*(Q[ei][idx_cell_wghost_y_BBB]  + Q[ei][idx_cell_wghost_y_TTT]) +
                             e_n*(Q[ei][idx_cell_wghost_y_BBBB] + Q[ei][idx_cell_wghost_y_TTTT])
                            )
                        );
                    }
                }
            }
            else if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*5)
            {
                for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_nghost = i + j*interior_dim_0;
                        
                        const int idx_cell_sound_speed = (i + num_subghosts_0_sound_speed) +
                            (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed;
                        
                        const int idx_cell_wghost_x_LLLLL = (i - 5 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_LLLL = (i - 4 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_LLL = (i - 3 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_LL = (i - 2 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_RR = (i + 2 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_RRR = (i + 3 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_RRRR = (i + 4 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_x_RRRRR = (i + 5 + num_subghosts_0_conservative_var) +
                            (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_BBBBB = (i + num_subghosts_0_conservative_var) +
                            (j - 5 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_BBBB = (i + num_subghosts_0_conservative_var) +
                            (j - 4 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_BBB = (i + num_subghosts_0_conservative_var) +
                            (j - 3 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_BB = (i + num_subghosts_0_conservative_var) +
                            (j - 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_B = (i + num_subghosts_0_conservative_var) +
                            (j - 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_T = (i + num_subghosts_0_conservative_var) +
                            (j + 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_TT = (i + num_subghosts_0_conservative_var) +
                            (j + 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_TTT = (i + num_subghosts_0_conservative_var) +
                            (j + 3 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_TTTT = (i + num_subghosts_0_conservative_var) +
                            (j + 4 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        const int idx_cell_wghost_y_TTTTT = (i + num_subghosts_0_conservative_var) +
                            (j + 5 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var;
                        
                        S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                            inv_dx_0*(
                             a_n*Q[ei][idx_cell_wghost] +
                             b_n*(Q[ei][idx_cell_wghost_x_L]     + Q[ei][idx_cell_wghost_x_R]) +
                             c_n*(Q[ei][idx_cell_wghost_x_LL]    + Q[ei][idx_cell_wghost_x_RR]) +
                             d_n*(Q[ei][idx_cell_wghost_x_LLL]   + Q[ei][idx_cell_wghost_x_RRR]) +
                             e_n*(Q[ei][idx_cell_wghost_x_LLLL]  + Q[ei][idx_cell_wghost_x_RRRR]) +
                             f_n*(Q[ei][idx_cell_wghost_x_LLLLL] + Q[ei][idx_cell_wghost_x_RRRRR])
                            ) +
                            inv_dx_1*(
                             a_n*Q[ei][idx_cell_wghost] +
                             b_n*(Q[ei][idx_cell_wghost_y_B]     + Q[ei][idx_cell_wghost_y_T]) +
                             c_n*(Q[ei][idx_cell_wghost_y_BB]    + Q[ei][idx_cell_wghost_y_TT]) +
                             d_n*(Q[ei][idx_cell_wghost_y_BBB]   + Q[ei][idx_cell_wghost_y_TTT]) +
                             e_n*(Q[ei][idx_cell_wghost_y_BBBB]  + Q[ei][idx_cell_wghost_y_TTTT]) +
                             f_n*(Q[ei][idx_cell_wghost_y_BBBBB] + Q[ei][idx_cell_wghost_y_TTTTT])
                            )
                        );
                    }
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the local lower indices and the number of cells in each dimension.
         */
        
        const int domain_lo_0 = domain_lo[0];
        const int domain_lo_1 = domain_lo[1];
        const int domain_lo_2 = domain_lo[2];
        const int domain_dim_0 = domain_dims[0];
        const int domain_dim_1 = domain_dims[1];
        const int domain_dim_2 = domain_dims[2];
        
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        /*
         * Get the grid spacings.
         */
        const Real inv_dx_0 = Real(1)/Real(dx[0]);
        const Real inv_dx_1 = Real(1)/Real(dx[1]);
        const Real inv_dx_2 = Real(1)/Real(dx[2]);
        
        /*
         * Get the number of ghost cells.
         */
        const int num_subghosts_0_sound_speed = num_subghosts_sound_speed[0];
        const int num_subghosts_1_sound_speed = num_subghosts_sound_speed[1];
        const int num_subghosts_2_sound_speed = num_subghosts_sound_speed[2];
        const int subghostcell_dim_0_sound_speed = subghostcell_dims_sound_speed[0];
        const int subghostcell_dim_1_sound_speed = subghostcell_dims_sound_speed[1];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Real* S = source->getPointer(ei);
            
            const int num_subghosts_0_conservative_var = num_subghosts_conservative_var[ei][0];
            const int num_subghosts_1_conservative_var = num_subghosts_conservative_var[ei][1];
            const int num_subghosts_2_conservative_var = num_subghosts_conservative_var[ei][2];
            const int subghostcell_dim_0_conservative_var = subghostcell_dims_conservative_var[ei][0];
            const int subghostcell_dim_1_conservative_var = subghostcell_dims_conservative_var[ei][1];
            
            if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*1)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_cell_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_B = (i + num_subghosts_0_conservative_var) +
                                (j - 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_T = (i + num_subghosts_0_conservative_var) +
                                (j + 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_B = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 1 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_F = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 1 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                                inv_dx_0*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_x_L] + Q[ei][idx_cell_wghost_x_R])
                                ) +
                                inv_dx_1*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_y_B] + Q[ei][idx_cell_wghost_y_T])
                                ) +
                                inv_dx_2*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_z_B] + Q[ei][idx_cell_wghost_z_F])
                                )
                            );
                        }
                    }
                }
            }
            else if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*2)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_cell_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_cell_wghost_x_LL = (i - 2 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_RR = (i + 2 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_BB = (i + num_subghosts_0_conservative_var) +
                                (j - 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_B = (i + num_subghosts_0_conservative_var) +
                                (j - 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_T = (i + num_subghosts_0_conservative_var) +
                                (j + 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_TT = (i + num_subghosts_0_conservative_var) +
                                (j + 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_BB = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 2 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_B = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 1 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_F = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 1 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_FF = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 2 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                                inv_dx_0*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_x_L]  + Q[ei][idx_cell_wghost_x_R]) +
                                 c_n*(Q[ei][idx_cell_wghost_x_LL] + Q[ei][idx_cell_wghost_x_RR])
                                ) +
                                inv_dx_1*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_y_B]  + Q[ei][idx_cell_wghost_y_T]) +
                                 c_n*(Q[ei][idx_cell_wghost_y_BB] + Q[ei][idx_cell_wghost_y_TT])
                                ) +
                                inv_dx_2*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_z_B]  + Q[ei][idx_cell_wghost_z_F]) +
                                 c_n*(Q[ei][idx_cell_wghost_z_BB] + Q[ei][idx_cell_wghost_z_FF])
                                )
                            );
                        }
                    }
                }
            }
            else if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*3)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_cell_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_cell_wghost_x_LLL = (i - 3 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_LL = (i - 2 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_RR = (i + 2 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_RRR = (i + 3 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_BBB = (i + num_subghosts_0_conservative_var) +
                                (j - 3 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_BB = (i + num_subghosts_0_conservative_var) +
                                (j - 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_B = (i + num_subghosts_0_conservative_var) +
                                (j - 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_T = (i + num_subghosts_0_conservative_var) +
                                (j + 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_TT = (i + num_subghosts_0_conservative_var) +
                                (j + 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_TTT = (i + num_subghosts_0_conservative_var) +
                                (j + 3 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_BBB = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 3 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_BB = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 2 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_B = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 1 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_F = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 1 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_FF = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 2 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_FFF = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 3 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                                inv_dx_0*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_x_L]   + Q[ei][idx_cell_wghost_x_R]) +
                                 c_n*(Q[ei][idx_cell_wghost_x_LL]  + Q[ei][idx_cell_wghost_x_RR]) +
                                 d_n*(Q[ei][idx_cell_wghost_x_LLL] + Q[ei][idx_cell_wghost_x_RRR])
                                ) +
                                inv_dx_1*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_y_B]   + Q[ei][idx_cell_wghost_y_T]) +
                                 c_n*(Q[ei][idx_cell_wghost_y_BB]  + Q[ei][idx_cell_wghost_y_TT]) +
                                 d_n*(Q[ei][idx_cell_wghost_y_BBB] + Q[ei][idx_cell_wghost_y_TTT])
                                ) +
                                inv_dx_2*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_z_B]   + Q[ei][idx_cell_wghost_z_F]) +
                                 c_n*(Q[ei][idx_cell_wghost_z_BB]  + Q[ei][idx_cell_wghost_z_FF]) +
                                 d_n*(Q[ei][idx_cell_wghost_z_BBB] + Q[ei][idx_cell_wghost_z_FFF])
                                )
                            );
                        }
                    }
                }
            }
            else if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*4)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_cell_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_cell_wghost_x_LLLL = (i - 4 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_LLL = (i - 3 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_LL = (i - 2 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_RR = (i + 2 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_RRR = (i + 3 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_RRRR = (i + 4 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_BBBB = (i + num_subghosts_0_conservative_var) +
                                (j - 4 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_BBB = (i + num_subghosts_0_conservative_var) +
                                (j - 3 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_BB = (i + num_subghosts_0_conservative_var) +
                                (j - 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_B = (i + num_subghosts_0_conservative_var) +
                                (j - 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_T = (i + num_subghosts_0_conservative_var) +
                                (j + 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_TT = (i + num_subghosts_0_conservative_var) +
                                (j + 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_TTT = (i + num_subghosts_0_conservative_var) +
                                (j + 3 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_TTTT = (i + num_subghosts_0_conservative_var) +
                                (j + 4 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_BBBB = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 4 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_BBB = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 3 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_BB = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 2 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_B = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 1 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_F = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 1 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_FF = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 2 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_FFF = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 3 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_FFFF = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 4 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                                inv_dx_0*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_x_L]    + Q[ei][idx_cell_wghost_x_R]) +
                                 c_n*(Q[ei][idx_cell_wghost_x_LL]   + Q[ei][idx_cell_wghost_x_RR]) +
                                 d_n*(Q[ei][idx_cell_wghost_x_LLL]  + Q[ei][idx_cell_wghost_x_RRR]) +
                                 e_n*(Q[ei][idx_cell_wghost_x_LLLL] + Q[ei][idx_cell_wghost_x_RRRR])
                                ) +
                                inv_dx_1*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_y_B]    + Q[ei][idx_cell_wghost_y_T]) +
                                 c_n*(Q[ei][idx_cell_wghost_y_BB]   + Q[ei][idx_cell_wghost_y_TT]) +
                                 d_n*(Q[ei][idx_cell_wghost_y_BBB]  + Q[ei][idx_cell_wghost_y_TTT]) +
                                 e_n*(Q[ei][idx_cell_wghost_y_BBBB] + Q[ei][idx_cell_wghost_y_TTTT])
                                ) +
                                inv_dx_2*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_z_B]    + Q[ei][idx_cell_wghost_z_F]) +
                                 c_n*(Q[ei][idx_cell_wghost_z_BB]   + Q[ei][idx_cell_wghost_z_FF]) +
                                 d_n*(Q[ei][idx_cell_wghost_z_BBB]  + Q[ei][idx_cell_wghost_z_FFF]) +
                                 e_n*(Q[ei][idx_cell_wghost_z_BBBB] + Q[ei][idx_cell_wghost_z_FFFF])
                                )
                            );
                        }
                    }
                }
            }
            else if (d_num_hyperviscosity_op_ghosts == hier::IntVector::getOne(d_dim)*5)
            {
                for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
                {
                    for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_cell_sound_speed = (i + num_subghosts_0_sound_speed) +
                                (j + num_subghosts_1_sound_speed)*subghostcell_dim_0_sound_speed +
                                (k + num_subghosts_2_sound_speed)*subghostcell_dim_0_sound_speed*
                                    subghostcell_dim_1_sound_speed;
                            
                            const int idx_cell_wghost_x_LLLLL = (i - 5 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_LLLL = (i - 4 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_LLL = (i - 3 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_LL = (i - 2 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_L = (i - 1 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_R = (i + 1 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_RR = (i + 2 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_RRR = (i + 3 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_RRRR = (i + 4 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_x_RRRRR = (i + 5 + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_BBBBB = (i + num_subghosts_0_conservative_var) +
                                (j - 5 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_BBBB = (i + num_subghosts_0_conservative_var) +
                                (j - 4 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_BBB = (i + num_subghosts_0_conservative_var) +
                                (j - 3 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_BB = (i + num_subghosts_0_conservative_var) +
                                (j - 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_B = (i + num_subghosts_0_conservative_var) +
                                (j - 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_T = (i + num_subghosts_0_conservative_var) +
                                (j + 1 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_TT = (i + num_subghosts_0_conservative_var) +
                                (j + 2 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_TTT = (i + num_subghosts_0_conservative_var) +
                                (j + 3 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_TTTT = (i + num_subghosts_0_conservative_var) +
                                (j + 4 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_y_TTTTT = (i + num_subghosts_0_conservative_var) +
                                (j + 5 + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_BBBBB = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 5 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_BBBB = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 4 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_BBB = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 3 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_BB = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 2 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_B = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k - 1 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_F = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 1 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_FF = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 2 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_FFF = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 3 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_FFFF = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 4 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            const int idx_cell_wghost_z_FFFFF = (i + num_subghosts_0_conservative_var) +
                                (j + num_subghosts_1_conservative_var)*subghostcell_dim_0_conservative_var +
                                (k + 5 + num_subghosts_2_conservative_var)*subghostcell_dim_0_conservative_var*
                                    subghostcell_dim_1_conservative_var;
                            
                            S[idx_cell_nghost] += Real(dt)*d_prefactor*c[idx_cell_sound_speed]*(
                                inv_dx_0*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_x_L]     + Q[ei][idx_cell_wghost_x_R]) +
                                 c_n*(Q[ei][idx_cell_wghost_x_LL]    + Q[ei][idx_cell_wghost_x_RR]) +
                                 d_n*(Q[ei][idx_cell_wghost_x_LLL]   + Q[ei][idx_cell_wghost_x_RRR]) +
                                 e_n*(Q[ei][idx_cell_wghost_x_LLLL]  + Q[ei][idx_cell_wghost_x_RRRR]) +
                                 f_n*(Q[ei][idx_cell_wghost_x_LLLLL] + Q[ei][idx_cell_wghost_x_RRRRR])
                                ) +
                                inv_dx_1*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_y_B]     + Q[ei][idx_cell_wghost_y_T]) +
                                 c_n*(Q[ei][idx_cell_wghost_y_BB]    + Q[ei][idx_cell_wghost_y_TT]) +
                                 d_n*(Q[ei][idx_cell_wghost_y_BBB]   + Q[ei][idx_cell_wghost_y_TTT]) +
                                 e_n*(Q[ei][idx_cell_wghost_y_BBBB]  + Q[ei][idx_cell_wghost_y_TTTT]) +
                                 f_n*(Q[ei][idx_cell_wghost_y_BBBBB] + Q[ei][idx_cell_wghost_y_TTTTT])
                                ) +
                                inv_dx_2*(
                                 a_n*Q[ei][idx_cell_wghost] +
                                 b_n*(Q[ei][idx_cell_wghost_z_B]     + Q[ei][idx_cell_wghost_z_F]) +
                                 c_n*(Q[ei][idx_cell_wghost_z_BB]    + Q[ei][idx_cell_wghost_z_FF]) +
                                 d_n*(Q[ei][idx_cell_wghost_z_BBB]   + Q[ei][idx_cell_wghost_z_FFF]) +
                                 e_n*(Q[ei][idx_cell_wghost_z_BBBB]  + Q[ei][idx_cell_wghost_z_FFFF]) +
                                 f_n*(Q[ei][idx_cell_wghost_z_BBBBB] + Q[ei][idx_cell_wghost_z_FFFFF])
                                )
                            );
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
}