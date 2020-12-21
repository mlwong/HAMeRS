#include "flow/convective_flux_reconstructors/DRP/ConvectiveFluxReconstructorDRP4.hpp"

/*
 * Timers interspersed throughout the class.
 */

HAMERS_SHARED_PTR<tbox::Timer> ConvectiveFluxReconstructorDRP4::t_reconstruct_flux;
HAMERS_SHARED_PTR<tbox::Timer> ConvectiveFluxReconstructorDRP4::t_compute_source;


ConvectiveFluxReconstructorDRP4::ConvectiveFluxReconstructorDRP4(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db):
        ConvectiveFluxReconstructor(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model,
            convective_flux_reconstructor_db)
{
    d_stencil_width = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("stencil_width", 13);
    d_stencil_width = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("d_stencil_width", d_stencil_width);
    
    if (d_stencil_width == 9)
    {
        d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*4;
    }
    else if (d_stencil_width == 11)
    {
        d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*5;
    }
    else if (d_stencil_width == 13)
    {
        d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*6;
    }
    else
    {
        TBOX_ERROR("ConvectiveFluxReconstructorDRP4::computeConvectiveFluxAndSourceOnPatch:"
            " Only 9-point, 11-point, 13-point stencil DRP schemes are implemented!");
    }
    
    d_eqn_form = d_flow_model->getEquationsForm();
    d_has_advective_eqn_form = false;
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
        {
            d_has_advective_eqn_form = true;
        }
    }
    
    t_reconstruct_flux = tbox::TimerManager::getManager()->
        getTimer("ConvectiveFluxReconstructorDRP4::t_reconstruct_flux");
    
    t_compute_source = tbox::TimerManager::getManager()->
        getTimer("ConvectiveFluxReconstructorDRP4::t_compute_source");
}


ConvectiveFluxReconstructorDRP4::~ConvectiveFluxReconstructorDRP4()
{
    t_reconstruct_flux.reset();
    t_compute_source.reset();
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorDRP4::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorDRP4 object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorDRP4: this = "
       << (ConvectiveFluxReconstructorDRP4 *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_stencil_width = "
       << d_stencil_width
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorDRP4::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putDouble("d_stencil_width", d_stencil_width);
}


/*
 * Compute the convective flux and source due to splitting of convective term on a patch.
 */
void
ConvectiveFluxReconstructorDRP4::computeConvectiveFluxAndSourceOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::SideVariable<double> >& variable_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_source,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(time);
    NULL_USE(RK_step_number);
    
    double a_n = double(0);
    double b_n = double(0);
    double c_n = double(0);
    double d_n = double(0);
    double e_n = double(0);
    double f_n = double(0);
    
    double a_m = double(0);
    double b_m = double(0);
    double c_m = double(0);
    double d_m = double(0);
    double e_m = double(0);
    double f_m = double(0);
    
    if (d_stencil_width = 9)
    {
        a_n = double( 0.841570125482);
        b_n = double(-0.244678631765);
        c_n = double( 0.059463584768);
        d_n = double(-0.007650904064);
        
        a_m = a_n + b_n + c_n + d_n;
        b_m = b_n + c_n + d_n;
        c_m = c_n + d_n;
        d_m = d_n;
    }
    else if (d_stencil_width = 11)
    {
        a_n = double( 0.872756993962);
        b_n = double(-0.286511173973);
        c_n = double( 0.090320001280);
        d_n = double(-0.020779405824);
        e_n = double( 0.002484594688);
        
        a_m = a_n + b_n + c_n + d_n + e_n;
        b_m = b_n + c_n + d_n + e_n;
        c_m = c_n + d_n + e_n;
        d_m = d_n + e_n;
        e_m = e_n;
    }
    else if (d_stencil_width = 13)
    {
        a_n = double( 0.907646591371);
        b_n = double(-0.337048393268);
        c_n = double( 0.133442885327);
        d_n = double(-0.045246480208);
        e_n = double( 0.011169294114);
        f_n = double(-0.001456501759);
        
        a_m = a_n + b_n + c_n + d_n + e_n + f_n;
        b_m = b_n + c_n + d_n + e_n + f_n;
        c_m = c_n + d_n + e_n + f_n;
        d_m = d_n + e_n + f_n;
        e_m = e_n + f_n;
        f_m = f_n;
    }
    else
    {
        TBOX_ERROR("ConvectiveFluxReconstructorDRP4::computeConvectiveFluxAndSourceOnPatch:"
            " Only 9-point, 11-point, 13-point stencil DRP schemes are implemented!");
    }
    
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
    HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux(
        HAMERS_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(variable_convective_flux, data_context)));
    
    // Get the cell data of source.
    HAMERS_SHARED_PTR<pdat::CellData<double> > source(
        HAMERS_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(variable_source, data_context)));
    
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    
    TBOX_ASSERT(source);
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // // Allocate temporary patch data.
    // HAMERS_SHARED_PTR<pdat::SideData<double> > velocity_midpoint;
    
    // if (d_has_advective_eqn_form)
    // {
    //     velocity_midpoint.reset(new pdat::SideData<double>(
    //         interior_box, d_dim.getValue(), hier::IntVector::getOne(d_dim)));
    // }
    
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
        
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        
        d_flow_model->registerDerivedVariables(num_subghosts_of_data);
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Get the pointers convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > convective_flux_node(2);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
        
        std::vector<double*> F_node_x;
        F_node_x.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
        }
        
       /*
         * Reconstruct the flux in the x-direction.
         */
        
        if (d_stencil_width = 9)
        {
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
                    
                    const int idx_node_LLLL = i - 4 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LLL  = i - 3 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LL   = i - 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_L    = i - 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_R    = i     + num_subghosts_0_convective_flux_x;
                    const int idx_node_RR   = i + 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRR  = i + 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRRR = i + 3 + num_subghosts_0_convective_flux_x;
                    
                    F_face_x[idx_face_x] = dt*(
                        a_m*(F_node_x[ei][idx_node_L]    + F_node_x[ei][idx_node_R]) +
                        b_m*(F_node_x[ei][idx_node_LL]   + F_node_x[ei][idx_node_RR]) +
                        c_m*(F_node_x[ei][idx_node_LLL]  + F_node_x[ei][idx_node_RRR]) +
                        d_m*(F_node_x[ei][idx_node_LLLL] + F_node_x[ei][idx_node_RRRR])
                        );
                }
            }
        }
        else if (d_stencil_width = 11)
        {
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
                    
                    const int idx_node_LLLLL = i - 5 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LLLL  = i - 4 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LLL   = i - 3 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LL    = i - 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_L     = i - 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_R     = i     + num_subghosts_0_convective_flux_x;
                    const int idx_node_RR    = i + 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRR   = i + 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRRR  = i + 3 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRRRR = i + 4 + num_subghosts_0_convective_flux_x;
                    
                    F_face_x[idx_face_x] = dt*(
                        a_m*(F_node_x[ei][idx_node_L]     + F_node_x[ei][idx_node_R]) +
                        b_m*(F_node_x[ei][idx_node_LL]    + F_node_x[ei][idx_node_RR]) +
                        c_m*(F_node_x[ei][idx_node_LLL]   + F_node_x[ei][idx_node_RRR]) +
                        d_m*(F_node_x[ei][idx_node_LLLL]  + F_node_x[ei][idx_node_RRRR]) +
                        e_m*(F_node_x[ei][idx_node_LLLLL] + F_node_x[ei][idx_node_RRRRR])
                        );
                }
            }
        }
        else if (d_stencil_width = 13)
        {
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
                    
                    const int idx_node_LLLLLL = i - 6 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LLLLL  = i - 5 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LLLL   = i - 4 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LLL    = i - 3 + num_subghosts_0_convective_flux_x;
                    const int idx_node_LL     = i - 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_L      = i - 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_R      = i     + num_subghosts_0_convective_flux_x;
                    const int idx_node_RR     = i + 1 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRR    = i + 2 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRRR   = i + 3 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRRRR  = i + 4 + num_subghosts_0_convective_flux_x;
                    const int idx_node_RRRRRR = i + 5 + num_subghosts_0_convective_flux_x;
                    
                    F_face_x[idx_face_x] = dt*(
                        a_m*(F_node_x[ei][idx_node_L]      + F_node_x[ei][idx_node_R]) +
                        b_m*(F_node_x[ei][idx_node_LL]     + F_node_x[ei][idx_node_RR]) +
                        c_m*(F_node_x[ei][idx_node_LLL]    + F_node_x[ei][idx_node_RRR]) +
                        d_m*(F_node_x[ei][idx_node_LLLL]   + F_node_x[ei][idx_node_RRRR]) +
                        e_m*(F_node_x[ei][idx_node_LLLLL]  + F_node_x[ei][idx_node_RRRRR]) +
                        f_m*(F_node_x[ei][idx_node_LLLLLL] + F_node_x[ei][idx_node_RRRRRR])
                        );
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
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
        
        d_flow_model->registerDerivedVariables(num_subghosts_of_data);
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Get the pointers convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > convective_flux_node(2);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
        
        const int num_subghosts_0_convective_flux_x = num_subghosts_convective_flux_x[0];
        const int num_subghosts_1_convective_flux_x = num_subghosts_convective_flux_x[1];
        const int subghostcell_dim_0_convective_flux_x = subghostcell_dims_convective_flux_x[0];
        
        const int num_subghosts_0_convective_flux_y = num_subghosts_convective_flux_y[0];
        const int num_subghosts_1_convective_flux_y = num_subghosts_convective_flux_y[1];
        const int subghostcell_dim_0_convective_flux_y = subghostcell_dims_convective_flux_y[0];
        
        std::vector<double*> F_node_x;
        std::vector<double*> F_node_y;
        F_node_x.reserve(d_num_eqn);
        F_node_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_node_x.push_back(convective_flux_node[0]->getPointer(ei));
            F_node_y.push_back(convective_flux_node[1]->getPointer(ei));
        }
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        if (d_stencil_width = 9)
        {
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
                        
                        const int idx_node_LLLL = (i - 4 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LLL  = (i - 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LL   = (i - 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_L    = (i - 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_R    = (i + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RR   = (i + 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRR  = (i + 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRRR = (i + 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        F_face_x[idx_face_x] = dt*(
                            a_m*(F_node_x[ei][idx_node_L]    + F_node_x[ei][idx_node_R]) +
                            b_m*(F_node_x[ei][idx_node_LL]   + F_node_x[ei][idx_node_RR]) +
                            c_m*(F_node_x[ei][idx_node_LLL]  + F_node_x[ei][idx_node_RRR]) +
                            d_m*(F_node_x[ei][idx_node_LLLL] + F_node_x[ei][idx_node_RRRR])
                            );
                    }
                }
            }
        }
        else if (d_stencil_width = 11)
        {
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
                        
                        const int idx_node_LLLLL = (i - 5 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LLLL  = (i - 4 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LLL   = (i - 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LL    = (i - 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_L     = (i - 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_R     = (i + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RR    = (i + 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRR   = (i + 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRRR  = (i + 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRRRR = (i + 4 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        F_face_x[idx_face_x] = dt*(
                            a_m*(F_node_x[ei][idx_node_L]     + F_node_x[ei][idx_node_R]) +
                            b_m*(F_node_x[ei][idx_node_LL]    + F_node_x[ei][idx_node_RR]) +
                            c_m*(F_node_x[ei][idx_node_LLL]   + F_node_x[ei][idx_node_RRR]) +
                            d_m*(F_node_x[ei][idx_node_LLLL]  + F_node_x[ei][idx_node_RRRR]) +
                            e_m*(F_node_x[ei][idx_node_LLLLL] + F_node_x[ei][idx_node_RRRRR])
                            );
                    }
                }
            }
        }
        else if (d_stencil_width = 13)
        {
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
                        
                        const int idx_node_LLLLLL = (i - 6 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LLLLL  = (i - 5 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LLLL   = (i - 4 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LLL    = (i - 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_LL     = (i - 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_L      = (i - 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_R      = (i + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RR     = (i + 1 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRR    = (i + 2 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRRR   = (i + 3 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRRRR  = (i + 4 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        const int idx_node_RRRRRR = (i + 5 + num_subghosts_0_convective_flux_x) +
                            (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x;
                        
                        F_face_x[idx_face_x] = dt*(
                            a_m*(F_node_x[ei][idx_node_L]      + F_node_x[ei][idx_node_R]) +
                            b_m*(F_node_x[ei][idx_node_LL]     + F_node_x[ei][idx_node_RR]) +
                            c_m*(F_node_x[ei][idx_node_LLL]    + F_node_x[ei][idx_node_RRR]) +
                            d_m*(F_node_x[ei][idx_node_LLLL]   + F_node_x[ei][idx_node_RRRR]) +
                            e_m*(F_node_x[ei][idx_node_LLLLL]  + F_node_x[ei][idx_node_RRRRR]) +
                            f_m*(F_node_x[ei][idx_node_LLLLLL] + F_node_x[ei][idx_node_RRRRRR])
                            );
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in the y-direction.
         */
        
        if (d_stencil_width = 9)
        {
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
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        const int idx_node_BBBB = (i + num_subghosts_0_convective_flux_y) +
                            (j - 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BBB  = (i + num_subghosts_0_convective_flux_y) +
                            (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BB   = (i + num_subghosts_0_convective_flux_y) +
                            (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_B    = (i + num_subghosts_0_convective_flux_y) +
                            (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_T    = (i + num_subghosts_0_convective_flux_y) +
                            (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TT   = (i + num_subghosts_0_convective_flux_y) +
                            (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTT  = (i + num_subghosts_0_convective_flux_y) +
                            (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTTT = (i + num_subghosts_0_convective_flux_y) +
                            (j + 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        F_face_y[idx_face_y] = dt*(
                            a_m*(F_node_y[ei][idx_node_B]    + F_node_y[ei][idx_node_T]) +
                            b_m*(F_node_y[ei][idx_node_BB]   + F_node_y[ei][idx_node_TT]) +
                            c_m*(F_node_y[ei][idx_node_BBB]  + F_node_y[ei][idx_node_TTT]) +
                            d_m*(F_node_y[ei][idx_node_BBBB] + F_node_y[ei][idx_node_TTTT])
                            );
                    }
                }
            }
        }
        else if (d_stencil_width = 11)
        {
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
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        const int idx_node_BBBBB = (i + num_subghosts_0_convective_flux_y) +
                            (j - 5 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BBBB  = (i + num_subghosts_0_convective_flux_y) +
                            (j - 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BBB   = (i + num_subghosts_0_convective_flux_y) +
                            (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BB    = (i + num_subghosts_0_convective_flux_y) +
                            (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_B     = (i + num_subghosts_0_convective_flux_y) +
                            (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_T     = (i + num_subghosts_0_convective_flux_y) +
                            (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TT    = (i + num_subghosts_0_convective_flux_y) +
                            (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTT   = (i + num_subghosts_0_convective_flux_y) +
                            (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTTT  = (i + num_subghosts_0_convective_flux_y) +
                            (j + 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTTTT = (i + num_subghosts_0_convective_flux_y) +
                            (j + 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        F_face_y[idx_face_y] = dt*(
                            a_m*(F_node_y[ei][idx_node_B]     + F_node_y[ei][idx_node_T]) +
                            b_m*(F_node_y[ei][idx_node_BB]    + F_node_y[ei][idx_node_TT]) +
                            c_m*(F_node_y[ei][idx_node_BBB]   + F_node_y[ei][idx_node_TTT]) +
                            d_m*(F_node_y[ei][idx_node_BBBB]  + F_node_y[ei][idx_node_TTTT]) +
                            e_m*(F_node_y[ei][idx_node_BBBBB] + F_node_y[ei][idx_node_TTTTT])
                            );
                    }
                }
            }
        }
        else if (d_stencil_width = 13)
        {
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
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        const int idx_node_BBBBBB = (i + num_subghosts_0_convective_flux_y) +
                            (j - 6 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BBBBB  = (i + num_subghosts_0_convective_flux_y) +
                            (j - 5 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BBBB   = (i + num_subghosts_0_convective_flux_y) +
                            (j - 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BBB    = (i + num_subghosts_0_convective_flux_y) +
                            (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_BB     = (i + num_subghosts_0_convective_flux_y) +
                            (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_B      = (i + num_subghosts_0_convective_flux_y) +
                            (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_T      = (i + num_subghosts_0_convective_flux_y) +
                            (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TT     = (i + num_subghosts_0_convective_flux_y) +
                            (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTT    = (i + num_subghosts_0_convective_flux_y) +
                            (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTTT   = (i + num_subghosts_0_convective_flux_y) +
                            (j + 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTTTT  = (i + num_subghosts_0_convective_flux_y) +
                            (j + 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        const int idx_node_TTTTTT = (i + num_subghosts_0_convective_flux_y) +
                            (j + 5 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y;
                        
                        F_face_y[idx_face_y] = dt*(
                            a_m*(F_node_y[ei][idx_node_B]      + F_node_y[ei][idx_node_T]) +
                            b_m*(F_node_y[ei][idx_node_BB]     + F_node_y[ei][idx_node_TT]) +
                            c_m*(F_node_y[ei][idx_node_BBB]    + F_node_y[ei][idx_node_TTT]) +
                            d_m*(F_node_y[ei][idx_node_BBBB]   + F_node_y[ei][idx_node_TTTT]) +
                            e_m*(F_node_y[ei][idx_node_BBBBB]  + F_node_y[ei][idx_node_TTTTT]) +
                            f_m*(F_node_y[ei][idx_node_BBBBBB] + F_node_y[ei][idx_node_TTTTTT])
                            );
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
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Z", d_num_conv_ghosts));
        
        d_flow_model->registerDerivedVariables(num_subghosts_of_data);
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Get the pointers convective flux cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > convective_flux_node(3);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
        convective_flux_node[2] = d_flow_model->getCellData("CONVECTIVE_FLUX_Z");
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_z = convective_flux_node[2]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_z = convective_flux_node[2]->getGhostBox().numberCells();
        
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
        
        /*
         * Reconstruct the flux in the x-direction.
         */
        
        if (d_stencil_width = 9)
        {
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
                                k*(interior_dim_0 + 1)*interior_dim_1;
                            
                            const int idx_node_LLLL = (i - 4 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LLL  = (i - 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LL   = (i - 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_L    = (i - 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_R    = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RR   = (i + 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRR  = (i + 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRRR = (i + 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            F_face_x[idx_face_x] = dt*(
                                a_m*(F_node_x[ei][idx_node_L]    + F_node_x[ei][idx_node_R]) +
                                b_m*(F_node_x[ei][idx_node_LL]   + F_node_x[ei][idx_node_RR]) +
                                c_m*(F_node_x[ei][idx_node_LLL]  + F_node_x[ei][idx_node_RRR]) +
                                d_m*(F_node_x[ei][idx_node_LLLL] + F_node_x[ei][idx_node_RRRR])
                                );
                        }
                    }
                }
            }
        }
        else if (d_stencil_width = 11)
        {
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
                                k*(interior_dim_0 + 1)*interior_dim_1;
                            
                            const int idx_node_LLLLL = (i - 5 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LLLL  = (i - 4 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LLL   = (i - 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LL    = (i - 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_L     = (i - 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_R     = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RR    = (i + 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRR   = (i + 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRRR  = (i + 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRRRR = (i + 4 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            F_face_x[idx_face_x] = dt*(
                                a_m*(F_node_x[ei][idx_node_L]     + F_node_x[ei][idx_node_R]) +
                                b_m*(F_node_x[ei][idx_node_LL]    + F_node_x[ei][idx_node_RR]) +
                                c_m*(F_node_x[ei][idx_node_LLL]   + F_node_x[ei][idx_node_RRR]) +
                                d_m*(F_node_x[ei][idx_node_LLLL]  + F_node_x[ei][idx_node_RRRR]) +
                                e_m*(F_node_x[ei][idx_node_LLLLL] + F_node_x[ei][idx_node_RRRRR])
                                );
                        }
                    }
                }
            }
        }
        else if (d_stencil_width = 13)
        {
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
                                k*(interior_dim_0 + 1)*interior_dim_1;
                            
                            const int idx_node_LLLLLL = (i - 6 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LLLLL  = (i - 5 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LLLL   = (i - 4 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LLL    = (i - 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_LL     = (i - 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_L      = (i - 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x + 
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_R      = (i + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RR     = (i + 1 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRR    = (i + 2 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRRR   = (i + 3 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRRRR  = (i + 4 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            const int idx_node_RRRRRR = (i + 5 + num_subghosts_0_convective_flux_x) +
                                (j + num_subghosts_1_convective_flux_x)*subghostcell_dim_0_convective_flux_x +
                                (k + num_subghosts_2_convective_flux_x)*subghostcell_dim_0_convective_flux_x*
                                    subghostcell_dim_1_convective_flux_x;
                            
                            F_face_x[idx_face_x] = dt*(
                                a_m*(F_node_x[ei][idx_node_L]      + F_node_x[ei][idx_node_R]) +
                                b_m*(F_node_x[ei][idx_node_LL]     + F_node_x[ei][idx_node_RR]) +
                                c_m*(F_node_x[ei][idx_node_LLL]    + F_node_x[ei][idx_node_RRR]) +
                                d_m*(F_node_x[ei][idx_node_LLLL]   + F_node_x[ei][idx_node_RRRR]) +
                                e_m*(F_node_x[ei][idx_node_LLLLL]  + F_node_x[ei][idx_node_RRRRR]) +
                                f_m*(F_node_x[ei][idx_node_LLLLLL] + F_node_x[ei][idx_node_RRRRRR])
                                );
                        }
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in the y-direction.
         */
        
        if (d_stencil_width = 9)
        {
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
                            const int idx_face_y = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            const int idx_node_BBBB = (i + num_subghosts_0_convective_flux_y) +
                                (j - 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BBB  = (i + num_subghosts_0_convective_flux_y) +
                                (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BB   = (i + num_subghosts_0_convective_flux_y) +
                                (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_B    = (i + num_subghosts_0_convective_flux_y) +
                                (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_T    = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TT   = (i + num_subghosts_0_convective_flux_y) +
                                (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTT  = (i + num_subghosts_0_convective_flux_y) +
                                (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTTT = (i + num_subghosts_0_convective_flux_y) +
                                (j + 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            F_face_y[idx_face_y] = dt*(
                                a_m*(F_node_y[ei][idx_node_B]    + F_node_y[ei][idx_node_T]) +
                                b_m*(F_node_y[ei][idx_node_BB]   + F_node_y[ei][idx_node_TT]) +
                                c_m*(F_node_y[ei][idx_node_BBB]  + F_node_y[ei][idx_node_TTT]) +
                                d_m*(F_node_y[ei][idx_node_BBBB] + F_node_y[ei][idx_node_TTTT])
                                );
                        }
                    }
                }
            }
        }
        else if (d_stencil_width = 11)
        {
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
                            const int idx_face_y = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            const int idx_node_BBBBB = (i + num_subghosts_0_convective_flux_y) +
                                (j - 5 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BBBB  = (i + num_subghosts_0_convective_flux_y) +
                                (j - 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BBB   = (i + num_subghosts_0_convective_flux_y) +
                                (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BB    = (i + num_subghosts_0_convective_flux_y) +
                                (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_B     = (i + num_subghosts_0_convective_flux_y) +
                                (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_T     = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TT    = (i + num_subghosts_0_convective_flux_y) +
                                (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTT   = (i + num_subghosts_0_convective_flux_y) +
                                (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTTT  = (i + num_subghosts_0_convective_flux_y) +
                                (j + 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTTTT = (i + num_subghosts_0_convective_flux_y) +
                                (j + 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            F_face_y[idx_face_y] = dt*(
                                a_m*(F_node_y[ei][idx_node_B]     + F_node_y[ei][idx_node_T]) +
                                b_m*(F_node_y[ei][idx_node_BB]    + F_node_y[ei][idx_node_TT]) +
                                c_m*(F_node_y[ei][idx_node_BBB]   + F_node_y[ei][idx_node_TTT]) +
                                d_m*(F_node_y[ei][idx_node_BBBB]  + F_node_y[ei][idx_node_TTTT]) +
                                e_m*(F_node_y[ei][idx_node_BBBBB] + F_node_y[ei][idx_node_TTTTT])
                                );
                        }
                    }
                }
            }
        }
        else if (d_stencil_width = 13)
        {
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
                            const int idx_face_y = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            const int idx_node_BBBBBB = (i + num_subghosts_0_convective_flux_y) +
                                (j - 6 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BBBBB  = (i + num_subghosts_0_convective_flux_y) +
                                (j - 5 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BBBB   = (i + num_subghosts_0_convective_flux_y) +
                                (j - 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BBB    = (i + num_subghosts_0_convective_flux_y) +
                                (j - 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_BB     = (i + num_subghosts_0_convective_flux_y) +
                                (j - 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_B      = (i + num_subghosts_0_convective_flux_y) +
                                (j - 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_T      = (i + num_subghosts_0_convective_flux_y) +
                                (j + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TT     = (i + num_subghosts_0_convective_flux_y) +
                                (j + 1 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTT    = (i + num_subghosts_0_convective_flux_y) +
                                (j + 2 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTTT   = (i + num_subghosts_0_convective_flux_y) +
                                (j + 3 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTTTT  = (i + num_subghosts_0_convective_flux_y) +
                                (j + 4 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            const int idx_node_TTTTTT = (i + num_subghosts_0_convective_flux_y) +
                                (j + 5 + num_subghosts_1_convective_flux_y)*subghostcell_dim_0_convective_flux_y +
                                (k + num_subghosts_2_convective_flux_y)*subghostcell_dim_0_convective_flux_y*
                                    subghostcell_dim_1_convective_flux_y;
                            
                            F_face_y[idx_face_y] = dt*(
                                a_m*(F_node_y[ei][idx_node_B]      + F_node_y[ei][idx_node_T]) +
                                b_m*(F_node_y[ei][idx_node_BB]     + F_node_y[ei][idx_node_TT]) +
                                c_m*(F_node_y[ei][idx_node_BBB]    + F_node_y[ei][idx_node_TTT]) +
                                d_m*(F_node_y[ei][idx_node_BBBB]   + F_node_y[ei][idx_node_TTTT]) +
                                e_m*(F_node_y[ei][idx_node_BBBBB]  + F_node_y[ei][idx_node_TTTTT]) +
                                f_m*(F_node_y[ei][idx_node_BBBBBB] + F_node_y[ei][idx_node_TTTTTT])
                                );
                        }
                    }
                }
            }
        }
        
        /*
         * Reconstruct the flux in the z-direction.
         */
        
        if (d_stencil_width = 9)
        {
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
                            const int idx_face_z = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_node_BBBB = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 4 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BBB  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BB   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_B    = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_F    = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FF   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFF  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFFF = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            F_face_z[idx_face_z] = dt*(
                                a_m*(F_node_z[ei][idx_node_B]    + F_node_z[ei][idx_node_F]) +
                                b_m*(F_node_z[ei][idx_node_BB]   + F_node_z[ei][idx_node_FF]) +
                                c_m*(F_node_z[ei][idx_node_BBB]  + F_node_z[ei][idx_node_FFF]) +
                                d_m*(F_node_z[ei][idx_node_BBBB] + F_node_z[ei][idx_node_FFFF])
                                );
                        }
                    }
                }
            }
        }
        else if (d_stencil_width = 11)
        {
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
                            const int idx_face_z = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_node_BBBBB = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 5 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BBBB  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 4 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BBB   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BB    = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_B     = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_F     = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FF    = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFF   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFFF  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFFFF = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 4 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            F_face_z[idx_face_z] = dt*(
                                a_m*(F_node_z[ei][idx_node_B]     + F_node_z[ei][idx_node_F]) +
                                b_m*(F_node_z[ei][idx_node_BB]    + F_node_z[ei][idx_node_FF]) +
                                c_m*(F_node_z[ei][idx_node_BBB]   + F_node_z[ei][idx_node_FFF]) +
                                d_m*(F_node_z[ei][idx_node_BBBB]  + F_node_z[ei][idx_node_FFFF]) +
                                e_m*(F_node_z[ei][idx_node_BBBBB] + F_node_z[ei][idx_node_FFFFF])
                                );
                        }
                    }
                }
            }
        }
        else if (d_stencil_width = 13)
        {
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
                            const int idx_face_z = i +
                                j*interior_dim_0 + 
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_node_BBBBBB = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 6 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BBBBB  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 5 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BBBB   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 4 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BBB    = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_BB     = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_B      = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k - 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_F      = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FF     = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 1 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFF    = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 2 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFFF   = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 3 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFFFF  = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 4 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            const int idx_node_FFFFFF = (i + num_subghosts_0_convective_flux_z) +
                                (j + num_subghosts_1_convective_flux_z)*subghostcell_dim_0_convective_flux_z +
                                (k + 5 + num_subghosts_2_convective_flux_z)*subghostcell_dim_0_convective_flux_z*
                                    subghostcell_dim_1_convective_flux_z;
                            
                            F_face_z[idx_face_z] = dt*(
                                a_m*(F_node_z[ei][idx_node_B]      + F_node_z[ei][idx_node_F]) +
                                b_m*(F_node_z[ei][idx_node_BB]     + F_node_z[ei][idx_node_FF]) +
                                c_m*(F_node_z[ei][idx_node_BBB]    + F_node_z[ei][idx_node_FFF]) +
                                d_m*(F_node_z[ei][idx_node_BBBB]   + F_node_z[ei][idx_node_FFFF]) +
                                e_m*(F_node_z[ei][idx_node_BBBBB]  + F_node_z[ei][idx_node_FFFFF]) +
                                f_m*(F_node_z[ei][idx_node_BBBBBB] + F_node_z[ei][idx_node_FFFFFF])
                                );
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
