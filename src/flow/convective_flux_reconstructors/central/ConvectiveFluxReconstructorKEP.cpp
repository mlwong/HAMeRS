#include "flow/convective_flux_reconstructors/central/ConvectiveFluxReconstructorKEP.hpp"

/*
 * Timers interspersed throughout the class.
 */

HAMERS_SHARED_PTR<tbox::Timer> ConvectiveFluxReconstructorKEP::t_reconstruct_flux;
HAMERS_SHARED_PTR<tbox::Timer> ConvectiveFluxReconstructorKEP::t_compute_source;


ConvectiveFluxReconstructorKEP::ConvectiveFluxReconstructorKEP(
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
    d_stencil_width = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("stencil_width", 3);
    d_stencil_width = d_convective_flux_reconstructor_db->
        getIntegerWithDefault("d_stencil_width", d_stencil_width);
    
    if (d_stencil_width == 3)
    {
        d_num_conv_ghosts = hier::IntVector::getOne(d_dim);
        
        d_coef_a = double(1)/double(2);
        d_coef_b = double(0);
        d_coef_c = double(0);
    }
    else if (d_stencil_width == 5)
    {
        d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*2;
        
        d_coef_a = double(2)/double(3);
        d_coef_b = -double(1)/double(12);
        d_coef_c = double(0);
    }
    else if (d_stencil_width == 7)
    {
        d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*3;
        
        d_coef_a = double(3)/double(4);
        d_coef_b = -double(3)/double(20);
        d_coef_c = double(1)/double(60);
    }
    else
    {
        TBOX_ERROR("ConvectiveFluxReconstructorKEP::computeConvectiveFluxAndSourceOnPatch:"
            " Only 3-point, 5-point, 7-point stencil central schemes are implemented!");
    }
    d_coef_d = double(0);
    d_coef_e = double(0);
    d_coef_f = double(0);
    
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
        getTimer("ConvectiveFluxReconstructorKEP::t_reconstruct_flux");
    
    t_compute_source = tbox::TimerManager::getManager()->
        getTimer("ConvectiveFluxReconstructorKEP::t_compute_source");
}


ConvectiveFluxReconstructorKEP::~ConvectiveFluxReconstructorKEP()
{
    t_reconstruct_flux.reset();
    t_compute_source.reset();
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorKEP::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorKEP object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorKEP: this = "
       << (ConvectiveFluxReconstructorKEP *)this
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
ConvectiveFluxReconstructorKEP::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putInteger("d_stencil_width", d_stencil_width);
}


/*
 * Compute the convective flux and source due to splitting of convective term on a patch.
 */
void
ConvectiveFluxReconstructorKEP::computeConvectiveFluxAndSourceOnPatch(
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
}


/*
 * Add linear term to convective flux in x-direction.
 */
void
ConvectiveFluxReconstructorKEP::addLinearTermToConvectiveFluxX(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_f,
    const int component_idx_flux,
    const int component_idx_f,
    const double dt) const
{
    const hier::Box interior_box = data_convective_flux->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_f->getBox().isSpatiallyEqual(interior_box));
#endif
    
    // Get the dimensions of the ghost cell box.
    const hier::Box ghost_box_f = data_f->getGhostBox();
    const hier::IntVector ghostcell_dims_f = ghost_box_f.numberCells();
    const hier::IntVector num_ghosts_f = data_f->getGhostCellWidth();
    
    /*
     * Get the pointers to the data.
     */
    
    double* F_face_x = data_convective_flux->getPointer(0, component_idx_flux);
    double* f        = data_f->getPointer(component_idx_f);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Get the number of ghost cells.
         */
        const int num_ghosts_0_f = num_ghosts_f[0];
        
        if (d_stencil_width == 3)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_L = i - 1 + num_ghosts_0_f;
                const int idx_f_R = i     + num_ghosts_0_f;
                
                F_face_x[idx_face_x] += dt*(
                    (d_coef_a*(f[idx_f_L] + f[idx_f_R])));
            }
        }
        else if (d_stencil_width == 5)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LL = i - 2 + num_ghosts_0_f;
                const int idx_f_L  = i - 1 + num_ghosts_0_f;
                const int idx_f_R  = i     + num_ghosts_0_f;
                const int idx_f_RR = i + 1 + num_ghosts_0_f;
                
                F_face_x[idx_face_x] += dt*(
                    (d_coef_a*(f[idx_f_L]  + f[idx_f_R] )) +
                    (d_coef_b*(f[idx_f_LL] + f[idx_f_R] )  +
                              (f[idx_f_L]  + f[idx_f_RR])));
            }
        }
        else if (d_stencil_width == 7)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LLL = i - 3 + num_ghosts_0_f;
                const int idx_f_LL  = i - 2 + num_ghosts_0_f;
                const int idx_f_L   = i - 1 + num_ghosts_0_f;
                const int idx_f_R   = i     + num_ghosts_0_f;
                const int idx_f_RR  = i + 1 + num_ghosts_0_f;
                const int idx_f_RRR = i + 2 + num_ghosts_0_f;
                
                F_face_x[idx_face_x] += dt*(
                    (d_coef_a*(f[idx_f_L]   + f[idx_f_R]  )) +
                    (d_coef_b*(f[idx_f_LL]  + f[idx_f_R]  )  +
                              (f[idx_f_L]   + f[idx_f_RR] )) +
                    (d_coef_c*(f[idx_f_LLL] + f[idx_f_R]  )  +
                              (f[idx_f_LL]  + f[idx_f_RR] )  +
                              (f[idx_f_L]   + f[idx_f_RRR])));
            }
        }
        else if (d_stencil_width == 9)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LLLL = i - 4 + num_ghosts_0_f;
                const int idx_f_LLL  = i - 3 + num_ghosts_0_f;
                const int idx_f_LL   = i - 2 + num_ghosts_0_f;
                const int idx_f_L    = i - 1 + num_ghosts_0_f;
                const int idx_f_R    = i     + num_ghosts_0_f;
                const int idx_f_RR   = i + 1 + num_ghosts_0_f;
                const int idx_f_RRR  = i + 2 + num_ghosts_0_f;
                const int idx_f_RRRR = i + 3 + num_ghosts_0_f;
                
                F_face_x[idx_face_x] += dt*(
                    (d_coef_a*(f[idx_f_L]    + f[idx_f_R]   )) +
                    (d_coef_b*(f[idx_f_LL]   + f[idx_f_R]   )  +
                              (f[idx_f_L]    + f[idx_f_RR]  )) +
                    (d_coef_c*(f[idx_f_LLL]  + f[idx_f_R]   )  +
                              (f[idx_f_LL]   + f[idx_f_RR]  )  +
                              (f[idx_f_L]    + f[idx_f_RRR] )) +
                    (d_coef_d*(f[idx_f_LLLL] + f[idx_f_R]   )  +
                              (f[idx_f_LLL]  + f[idx_f_RR]  )  +
                              (f[idx_f_LL]   + f[idx_f_RRR] )  +
                              (f[idx_f_L]    + f[idx_f_RRRR])));
            }
        }
        else if (d_stencil_width == 11)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LLLLL = i - 5 + num_ghosts_0_f;
                const int idx_f_LLLL  = i - 4 + num_ghosts_0_f;
                const int idx_f_LLL   = i - 3 + num_ghosts_0_f;
                const int idx_f_LL    = i - 2 + num_ghosts_0_f;
                const int idx_f_L     = i - 1 + num_ghosts_0_f;
                const int idx_f_R     = i     + num_ghosts_0_f;
                const int idx_f_RR    = i + 1 + num_ghosts_0_f;
                const int idx_f_RRR   = i + 2 + num_ghosts_0_f;
                const int idx_f_RRRR  = i + 3 + num_ghosts_0_f;
                const int idx_f_RRRRR = i + 4 + num_ghosts_0_f;
                
                F_face_x[idx_face_x] += dt*(
                    (d_coef_a*(f[idx_f_L]     + f[idx_f_R]    )) +
                    (d_coef_b*(f[idx_f_LL]    + f[idx_f_R]    )  +
                              (f[idx_f_L]     + f[idx_f_RR]   )) +
                    (d_coef_c*(f[idx_f_LLL]   + f[idx_f_R]    )  +
                              (f[idx_f_LL]    + f[idx_f_RR]   )  +
                              (f[idx_f_L]     + f[idx_f_RRR]  )) +
                    (d_coef_d*(f[idx_f_LLLL]  + f[idx_f_R]    )  +
                              (f[idx_f_LLL]   + f[idx_f_RR]   )  +
                              (f[idx_f_LL]    + f[idx_f_RRR]  )  +
                              (f[idx_f_L]     + f[idx_f_RRRR] )) +
                    (d_coef_e*(f[idx_f_LLLLL] + f[idx_f_R]    )  +
                              (f[idx_f_LLLL]  + f[idx_f_RR]   )  +
                              (f[idx_f_LLL]   + f[idx_f_RRR]  )  +
                              (f[idx_f_LL]    + f[idx_f_RRRR] )  +
                              (f[idx_f_L]     + f[idx_f_RRRRR])));
            }
        }
        else if (d_stencil_width == 13)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LLLLLL = i - 6 + num_ghosts_0_f;
                const int idx_f_LLLLL  = i - 5 + num_ghosts_0_f;
                const int idx_f_LLLL   = i - 4 + num_ghosts_0_f;
                const int idx_f_LLL    = i - 3 + num_ghosts_0_f;
                const int idx_f_LL     = i - 2 + num_ghosts_0_f;
                const int idx_f_L      = i - 1 + num_ghosts_0_f;
                const int idx_f_R      = i     + num_ghosts_0_f;
                const int idx_f_RR     = i + 1 + num_ghosts_0_f;
                const int idx_f_RRR    = i + 2 + num_ghosts_0_f;
                const int idx_f_RRRR   = i + 3 + num_ghosts_0_f;
                const int idx_f_RRRRR  = i + 4 + num_ghosts_0_f;
                const int idx_f_RRRRRR = i + 5 + num_ghosts_0_f;
                
                F_face_x[idx_face_x] += dt*(
                    (d_coef_a*(f[idx_f_L]      + f[idx_f_R]     )) +
                    (d_coef_b*(f[idx_f_LL]     + f[idx_f_R]     )  +
                              (f[idx_f_L]      + f[idx_f_RR]    )) +
                    (d_coef_c*(f[idx_f_LLL]    + f[idx_f_R]     )  +
                              (f[idx_f_LL]     + f[idx_f_RR]    )  +
                              (f[idx_f_L]      + f[idx_f_RRR]   )) +
                    (d_coef_d*(f[idx_f_LLLL]   + f[idx_f_R]     )  +
                              (f[idx_f_LLL]    + f[idx_f_RR]    )  +
                              (f[idx_f_LL]     + f[idx_f_RRR]   )  +
                              (f[idx_f_L]      + f[idx_f_RRRR]  )) +
                    (d_coef_e*(f[idx_f_LLLLL]  + f[idx_f_R]     )  +
                              (f[idx_f_LLLL]   + f[idx_f_RR]    )  +
                              (f[idx_f_LLL]    + f[idx_f_RRR]   )  +
                              (f[idx_f_LL]     + f[idx_f_RRRR]  )  +
                              (f[idx_f_L]      + f[idx_f_RRRRR] )) +
                    (d_coef_f*(f[idx_f_LLLLLL] + f[idx_f_R]     )  +
                              (f[idx_f_LLLLL]  + f[idx_f_RR]    )  +
                              (f[idx_f_LLLL]   + f[idx_f_RRR]   )  +
                              (f[idx_f_LLL]    + f[idx_f_RRRR]  )  +
                              (f[idx_f_LL]     + f[idx_f_RRRRR] )  +
                              (f[idx_f_L]      + f[idx_f_RRRRRR])));
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        /*
         * Get the number of ghost cells and the dimensions of the ghost cell box.
         */
        
        const int num_ghosts_0_f = num_ghosts_f[0];
        const int num_ghosts_1_f = num_ghosts_f[1];
        const int ghostcell_dim_0_f = ghostcell_dims_f[0];
        
        if (d_stencil_width == 3)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_L = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    F_face_x[idx_face_x] += dt*(
                        (d_coef_a*(f[idx_f_L] + f[idx_f_R])));
                }
            }
        }
        else if (d_stencil_width == 5)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LL = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L  = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R  = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    F_face_x[idx_face_x] += dt*(
                        (d_coef_a*(f[idx_f_L]  + f[idx_f_R] )) +
                        (d_coef_b*(f[idx_f_LL] + f[idx_f_R] )  +
                                  (f[idx_f_L]  + f[idx_f_RR])));
                }
            }
        }
        else if (d_stencil_width == 7)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LLL = (i - 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LL  = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L   = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R   = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR  = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRR = (i + 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    F_face_x[idx_face_x] += dt*(
                        (d_coef_a*(f[idx_f_L]   + f[idx_f_R]  )) +
                        (d_coef_b*(f[idx_f_LL]  + f[idx_f_R]  )  +
                                  (f[idx_f_L]   + f[idx_f_RR] )) +
                        (d_coef_c*(f[idx_f_LLL] + f[idx_f_R]  )  +
                                  (f[idx_f_LL]  + f[idx_f_RR] )  +
                                  (f[idx_f_L]   + f[idx_f_RRR])));
                }
            }
        }
        else if (d_stencil_width == 9)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LLLL = (i - 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLL  = (i - 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LL   = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L    = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R    = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR   = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRR  = (i + 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRR = (i + 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    F_face_x[idx_face_x] += dt*(
                        (d_coef_a*(f[idx_f_L]    + f[idx_f_R]   )) +
                        (d_coef_b*(f[idx_f_LL]   + f[idx_f_R]   )  +
                                  (f[idx_f_L]    + f[idx_f_RR]  )) +
                        (d_coef_c*(f[idx_f_LLL]  + f[idx_f_R]   )  +
                                  (f[idx_f_LL]   + f[idx_f_RR]  )  +
                                  (f[idx_f_L]    + f[idx_f_RRR] )) +
                        (d_coef_d*(f[idx_f_LLLL] + f[idx_f_R]   )  +
                                  (f[idx_f_LLL]  + f[idx_f_RR]  )  +
                                  (f[idx_f_LL]   + f[idx_f_RRR] )  +
                                  (f[idx_f_L]    + f[idx_f_RRRR])));
                }
            }
        }
        else if (d_stencil_width == 11)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LLLLL = (i - 5 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLLL  = (i - 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLL   = (i - 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LL    = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L     = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R     = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR    = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRR   = (i + 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRR  = (i + 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRRR = (i + 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    F_face_x[idx_face_x] += dt*(
                        (d_coef_a*(f[idx_f_L]     + f[idx_f_R]    )) +
                        (d_coef_b*(f[idx_f_LL]    + f[idx_f_R]    )  +
                                  (f[idx_f_L]     + f[idx_f_RR]   )) +
                        (d_coef_c*(f[idx_f_LLL]   + f[idx_f_R]    )  +
                                  (f[idx_f_LL]    + f[idx_f_RR]   )  +
                                  (f[idx_f_L]     + f[idx_f_RRR]  )) +
                        (d_coef_d*(f[idx_f_LLLL]  + f[idx_f_R]    )  +
                                  (f[idx_f_LLL]   + f[idx_f_RR]   )  +
                                  (f[idx_f_LL]    + f[idx_f_RRR]  )  +
                                  (f[idx_f_L]     + f[idx_f_RRRR] )) +
                        (d_coef_e*(f[idx_f_LLLLL] + f[idx_f_R]    )  +
                                  (f[idx_f_LLLL]  + f[idx_f_RR]   )  +
                                  (f[idx_f_LLL]   + f[idx_f_RRR]  )  +
                                  (f[idx_f_LL]    + f[idx_f_RRRR] )  +
                                  (f[idx_f_L]     + f[idx_f_RRRRR])));
                }
            }
        }
        else if (d_stencil_width == 13)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LLLLLL = (i - 6 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLLLL  = (i - 5 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLLL   = (i - 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLL    = (i - 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LL     = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L      = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R      = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR     = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRR    = (i + 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRR   = (i + 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRRR  = (i + 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRRRR = (i + 5 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    F_face_x[idx_face_x] += dt*(
                        (d_coef_a*(f[idx_f_L]      + f[idx_f_R]     )) +
                        (d_coef_b*(f[idx_f_LL]     + f[idx_f_R]     )  +
                                  (f[idx_f_L]      + f[idx_f_RR]    )) +
                        (d_coef_c*(f[idx_f_LLL]    + f[idx_f_R]     )  +
                                  (f[idx_f_LL]     + f[idx_f_RR]    )  +
                                  (f[idx_f_L]      + f[idx_f_RRR]   )) +
                        (d_coef_d*(f[idx_f_LLLL]   + f[idx_f_R]     )  +
                                  (f[idx_f_LLL]    + f[idx_f_RR]    )  +
                                  (f[idx_f_LL]     + f[idx_f_RRR]   )  +
                                  (f[idx_f_L]      + f[idx_f_RRRR]  )) +
                        (d_coef_e*(f[idx_f_LLLLL]  + f[idx_f_R]     )  +
                                  (f[idx_f_LLLL]   + f[idx_f_RR]    )  +
                                  (f[idx_f_LLL]    + f[idx_f_RRR]   )  +
                                  (f[idx_f_LL]     + f[idx_f_RRRR]  )  +
                                  (f[idx_f_L]      + f[idx_f_RRRRR] )) +
                        (d_coef_f*(f[idx_f_LLLLLL] + f[idx_f_R]     )  +
                                  (f[idx_f_LLLLL]  + f[idx_f_RR]    )  +
                                  (f[idx_f_LLLL]   + f[idx_f_RRR]   )  +
                                  (f[idx_f_LLL]    + f[idx_f_RRRR]  )  +
                                  (f[idx_f_LL]     + f[idx_f_RRRRR] )  +
                                  (f[idx_f_L]      + f[idx_f_RRRRRR])));
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        /*
         * Get the number of ghost cells and the dimensions of the ghost cell box.
         */
        
        const int num_ghosts_0_f = num_ghosts_f[0];
        const int num_ghosts_1_f = num_ghosts_f[1];
        const int num_ghosts_2_f = num_ghosts_f[2];
        const int ghostcell_dim_0_f = ghostcell_dims_f[0];
        const int ghostcell_dim_1_f = ghostcell_dims_f[1];
        
        if (d_stencil_width == 3)
        {
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
                        
                        const int idx_f_L = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_x[idx_face_x] += dt*(
                            (d_coef_a*(f[idx_f_L] + f[idx_f_R])));
                    }
                }
            }
        }
        else if (d_stencil_width == 5)
        {
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
                        
                        const int idx_f_LL = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L  = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_x[idx_face_x] += dt*(
                            (d_coef_a*(f[idx_f_L]  + f[idx_f_R] )) +
                            (d_coef_b*(f[idx_f_LL] + f[idx_f_R] )  +
                                      (f[idx_f_L]  + f[idx_f_RR])));
                    }
                }
            }
        }
        else if (d_stencil_width == 7)
        {
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
                        
                        const int idx_f_LLL = (i - 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LL  = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L   = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR  = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRR = (i + 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_x[idx_face_x] += dt*(
                            (d_coef_a*(f[idx_f_L]   + f[idx_f_R]  )) +
                            (d_coef_b*(f[idx_f_LL]  + f[idx_f_R]  )  +
                                      (f[idx_f_L]   + f[idx_f_RR] )) +
                            (d_coef_c*(f[idx_f_LLL] + f[idx_f_R]  )  +
                                      (f[idx_f_LL]  + f[idx_f_RR] )  +
                                      (f[idx_f_L]   + f[idx_f_RRR])));
                    }
                }
            }
        }
        else if (d_stencil_width == 9)
        {
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
                        
                        const int idx_f_LLLL = (i - 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLL  = (i - 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LL   = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L    = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR   = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRR  = (i + 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRR = (i + 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_x[idx_face_x] += dt*(
                            (d_coef_a*(f[idx_f_L]    + f[idx_f_R]   )) +
                            (d_coef_b*(f[idx_f_LL]   + f[idx_f_R]   )  +
                                      (f[idx_f_L]    + f[idx_f_RR]  )) +
                            (d_coef_c*(f[idx_f_LLL]  + f[idx_f_R]   )  +
                                      (f[idx_f_LL]   + f[idx_f_RR]  )  +
                                      (f[idx_f_L]    + f[idx_f_RRR] )) +
                            (d_coef_d*(f[idx_f_LLLL] + f[idx_f_R]   )  +
                                      (f[idx_f_LLL]  + f[idx_f_RR]  )  +
                                      (f[idx_f_LL]   + f[idx_f_RRR] )  +
                                      (f[idx_f_L]    + f[idx_f_RRRR])));
                    }
                }
            }
        }
        else if (d_stencil_width == 11)
        {
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
                        
                        const int idx_f_LLLLL = (i - 5 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLLL  = (i - 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLL   = (i - 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LL    = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L     = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR    = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRR   = (i + 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRR  = (i + 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRRR = (i + 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_x[idx_face_x] += dt*(
                            (d_coef_a*(f[idx_f_L]     + f[idx_f_R]    )) +
                            (d_coef_b*(f[idx_f_LL]    + f[idx_f_R]    )  +
                                      (f[idx_f_L]     + f[idx_f_RR]   )) +
                            (d_coef_c*(f[idx_f_LLL]   + f[idx_f_R]    )  +
                                      (f[idx_f_LL]    + f[idx_f_RR]   )  +
                                      (f[idx_f_L]     + f[idx_f_RRR]  )) +
                            (d_coef_d*(f[idx_f_LLLL]  + f[idx_f_R]    )  +
                                      (f[idx_f_LLL]   + f[idx_f_RR]   )  +
                                      (f[idx_f_LL]    + f[idx_f_RRR]  )  +
                                      (f[idx_f_L]     + f[idx_f_RRRR] )) +
                            (d_coef_e*(f[idx_f_LLLLL] + f[idx_f_R]    )  +
                                      (f[idx_f_LLLL]  + f[idx_f_RR]   )  +
                                      (f[idx_f_LLL]   + f[idx_f_RRR]  )  +
                                      (f[idx_f_LL]    + f[idx_f_RRRR] )  +
                                      (f[idx_f_L]     + f[idx_f_RRRRR])));
                    }
                }
            }
        }
        else if (d_stencil_width == 13)
        {
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
                        
                        const int idx_f_LLLLLL = (i - 6 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLLLL  = (i - 5 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLLL   = (i - 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLL    = (i - 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LL     = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L      = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R      = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR     = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRR    = (i + 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRR   = (i + 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRRR  = (i + 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRRRR = (i + 5 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_x[idx_face_x] += dt*(
                            (d_coef_a*(f[idx_f_L]      + f[idx_f_R]     )) +
                            (d_coef_b*(f[idx_f_LL]     + f[idx_f_R]     )  +
                                      (f[idx_f_L]      + f[idx_f_RR]    )) +
                            (d_coef_c*(f[idx_f_LLL]    + f[idx_f_R]     )  +
                                      (f[idx_f_LL]     + f[idx_f_RR]    )  +
                                      (f[idx_f_L]      + f[idx_f_RRR]   )) +
                            (d_coef_d*(f[idx_f_LLLL]   + f[idx_f_R]     )  +
                                      (f[idx_f_LLL]    + f[idx_f_RR]    )  +
                                      (f[idx_f_LL]     + f[idx_f_RRR]   )  +
                                      (f[idx_f_L]      + f[idx_f_RRRR]  )) +
                            (d_coef_e*(f[idx_f_LLLLL]  + f[idx_f_R]     )  +
                                      (f[idx_f_LLLL]   + f[idx_f_RR]    )  +
                                      (f[idx_f_LLL]    + f[idx_f_RRR]   )  +
                                      (f[idx_f_LL]     + f[idx_f_RRRR]  )  +
                                      (f[idx_f_L]      + f[idx_f_RRRRR] )) +
                            (d_coef_f*(f[idx_f_LLLLLL] + f[idx_f_R]     )  +
                                      (f[idx_f_LLLLL]  + f[idx_f_RR]    )  +
                                      (f[idx_f_LLLL]   + f[idx_f_RRR]   )  +
                                      (f[idx_f_LLL]    + f[idx_f_RRRR]  )  +
                                      (f[idx_f_LL]     + f[idx_f_RRRRR] )  +
                                      (f[idx_f_L]      + f[idx_f_RRRRRR])));
                    }
                }
            }
        }
    }
}


/*
 * Add quadratic term to convective flux in x-direction.
 */
void
ConvectiveFluxReconstructorKEP::addQuadraticTermToConvectiveFluxX(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_f,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_g,
    const int component_idx_flux,
    const int component_idx_f,
    const int component_idx_g,
    const double dt) const
{
    const hier::Box interior_box = data_convective_flux->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_f->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_g->getBox().isSpatiallyEqual(interior_box));
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_f = data_f->getGhostBox();
    const hier::IntVector ghostcell_dims_f = ghost_box_f.numberCells();
    const hier::IntVector num_ghosts_f = data_f->getGhostCellWidth();
    
    const hier::Box ghost_box_g = data_g->getGhostBox();
    const hier::IntVector ghostcell_dims_g = ghost_box_g.numberCells();
    const hier::IntVector num_ghosts_g = data_g->getGhostCellWidth();
    
    /*
     * Get the pointers to the data.
     */
    
    double* F_face_x = data_convective_flux->getPointer(0, component_idx_flux);
    double* f        = data_f->getPointer(component_idx_f);
    double* g        = data_g->getPointer(component_idx_g);
    
    const double half = double(1)/double(2);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Get the numbers of ghost cells.
         */
        const int num_ghosts_0_f = num_ghosts_f[0];
        const int num_ghosts_0_g = num_ghosts_g[0];
        
        if (d_stencil_width == 3)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_L = i - 1 + num_ghosts_0_f;
                const int idx_f_R = i     + num_ghosts_0_f;
                
                const int idx_g_L = i - 1 + num_ghosts_0_g;
                const int idx_g_R = i     + num_ghosts_0_g;
                
                F_face_x[idx_face_x] += dt*(
                    half*(d_coef_a*(f[idx_f_L] + f[idx_f_R])*(g[idx_g_L] + g[idx_g_R])));
            }
        }
        else if (d_stencil_width == 5)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LL = i - 2 + num_ghosts_0_f;
                const int idx_f_L  = i - 1 + num_ghosts_0_f;
                const int idx_f_R  = i     + num_ghosts_0_f;
                const int idx_f_RR = i + 1 + num_ghosts_0_f;
                
                const int idx_g_LL = i - 2 + num_ghosts_0_g;
                const int idx_g_L  = i - 1 + num_ghosts_0_g;
                const int idx_g_R  = i     + num_ghosts_0_g;
                const int idx_g_RR = i + 1 + num_ghosts_0_g;
                
                F_face_x[idx_face_x] += dt*(
                    half*(d_coef_a*(f[idx_f_L]  + f[idx_f_R] )*(g[idx_g_L]  + g[idx_g_R] )) +
                    half*(d_coef_b*(f[idx_f_LL] + f[idx_f_R] )*(g[idx_g_LL] + g[idx_g_R] )  +
                                   (f[idx_f_L]  + f[idx_f_RR])*(g[idx_g_L]  + g[idx_g_RR])));
            }
        }
        else if (d_stencil_width == 7)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LLL = i - 3 + num_ghosts_0_f;
                const int idx_f_LL  = i - 2 + num_ghosts_0_f;
                const int idx_f_L   = i - 1 + num_ghosts_0_f;
                const int idx_f_R   = i     + num_ghosts_0_f;
                const int idx_f_RR  = i + 1 + num_ghosts_0_f;
                const int idx_f_RRR = i + 2 + num_ghosts_0_f;
                
                const int idx_g_LLL = i - 3 + num_ghosts_0_g;
                const int idx_g_LL  = i - 2 + num_ghosts_0_g;
                const int idx_g_L   = i - 1 + num_ghosts_0_g;
                const int idx_g_R   = i     + num_ghosts_0_g;
                const int idx_g_RR  = i + 1 + num_ghosts_0_g;
                const int idx_g_RRR = i + 2 + num_ghosts_0_g;
                
                F_face_x[idx_face_x] += dt*(
                    half*(d_coef_a*(f[idx_f_L]   + f[idx_f_R]  )*(g[idx_g_L]   + g[idx_g_R]  )) +
                    half*(d_coef_b*(f[idx_f_LL]  + f[idx_f_R]  )*(g[idx_g_LL]  + g[idx_g_R]  )  +
                                   (f[idx_f_L]   + f[idx_f_RR] )*(g[idx_g_L]   + g[idx_g_RR] )) +
                    half*(d_coef_c*(f[idx_f_LLL] + f[idx_f_R]  )*(g[idx_g_LLL] + g[idx_g_R]  )  +
                                   (f[idx_f_LL]  + f[idx_f_RR] )*(g[idx_g_LL]  + g[idx_g_RR] )  +
                                   (f[idx_f_L]   + f[idx_f_RRR])*(g[idx_g_L]   + g[idx_g_RRR])));
            }
        }
        else if (d_stencil_width == 9)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LLLL = i - 4 + num_ghosts_0_f;
                const int idx_f_LLL  = i - 3 + num_ghosts_0_f;
                const int idx_f_LL   = i - 2 + num_ghosts_0_f;
                const int idx_f_L    = i - 1 + num_ghosts_0_f;
                const int idx_f_R    = i     + num_ghosts_0_f;
                const int idx_f_RR   = i + 1 + num_ghosts_0_f;
                const int idx_f_RRR  = i + 2 + num_ghosts_0_f;
                const int idx_f_RRRR = i + 3 + num_ghosts_0_f;
                
                const int idx_g_LLLL = i - 4 + num_ghosts_0_g;
                const int idx_g_LLL  = i - 3 + num_ghosts_0_g;
                const int idx_g_LL   = i - 2 + num_ghosts_0_g;
                const int idx_g_L    = i - 1 + num_ghosts_0_g;
                const int idx_g_R    = i     + num_ghosts_0_g;
                const int idx_g_RR   = i + 1 + num_ghosts_0_g;
                const int idx_g_RRR  = i + 2 + num_ghosts_0_g;
                const int idx_g_RRRR = i + 3 + num_ghosts_0_g;
                
                F_face_x[idx_face_x] += dt*(
                    half*(d_coef_a*(f[idx_f_L]    + f[idx_f_R]   )*(g[idx_g_L]    + g[idx_g_R]   )) +
                    half*(d_coef_b*(f[idx_f_LL]   + f[idx_f_R]   )*(g[idx_g_LL]   + g[idx_g_R]   )  +
                                   (f[idx_f_L]    + f[idx_f_RR]  )*(g[idx_g_L]    + g[idx_g_RR]  )) +
                    half*(d_coef_c*(f[idx_f_LLL]  + f[idx_f_R]   )*(g[idx_g_LLL]  + g[idx_g_R]   )  +
                                   (f[idx_f_LL]   + f[idx_f_RR]  )*(g[idx_g_LL]   + g[idx_g_RR]  )  +
                                   (f[idx_f_L]    + f[idx_f_RRR] )*(g[idx_g_L]    + g[idx_g_RRR] )) +
                    half*(d_coef_d*(f[idx_f_LLLL] + f[idx_f_R]   )*(g[idx_g_LLLL] + g[idx_g_R]   )  +
                                   (f[idx_f_LLL]  + f[idx_f_RR]  )*(g[idx_g_LLL]  + g[idx_g_RR]  )  +
                                   (f[idx_f_LL]   + f[idx_f_RRR] )*(g[idx_g_LL]   + g[idx_g_RRR] )  +
                                   (f[idx_f_L]    + f[idx_f_RRRR])*(g[idx_g_L]    + g[idx_g_RRRR])));
            }
        }
        else if (d_stencil_width == 11)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LLLLL = i - 5 + num_ghosts_0_f;
                const int idx_f_LLLL  = i - 4 + num_ghosts_0_f;
                const int idx_f_LLL   = i - 3 + num_ghosts_0_f;
                const int idx_f_LL    = i - 2 + num_ghosts_0_f;
                const int idx_f_L     = i - 1 + num_ghosts_0_f;
                const int idx_f_R     = i     + num_ghosts_0_f;
                const int idx_f_RR    = i + 1 + num_ghosts_0_f;
                const int idx_f_RRR   = i + 2 + num_ghosts_0_f;
                const int idx_f_RRRR  = i + 3 + num_ghosts_0_f;
                const int idx_f_RRRRR = i + 4 + num_ghosts_0_f;
                
                const int idx_g_LLLLL = i - 5 + num_ghosts_0_g;
                const int idx_g_LLLL  = i - 4 + num_ghosts_0_g;
                const int idx_g_LLL   = i - 3 + num_ghosts_0_g;
                const int idx_g_LL    = i - 2 + num_ghosts_0_g;
                const int idx_g_L     = i - 1 + num_ghosts_0_g;
                const int idx_g_R     = i     + num_ghosts_0_g;
                const int idx_g_RR    = i + 1 + num_ghosts_0_g;
                const int idx_g_RRR   = i + 2 + num_ghosts_0_g;
                const int idx_g_RRRR  = i + 3 + num_ghosts_0_g;
                const int idx_g_RRRRR = i + 4 + num_ghosts_0_g;
                
                F_face_x[idx_face_x] += dt*(
                    half*(d_coef_a*(f[idx_f_L]     + f[idx_f_R]    )*(g[idx_g_L]     + g[idx_g_R]    )) +
                    half*(d_coef_b*(f[idx_f_LL]    + f[idx_f_R]    )*(g[idx_g_LL]    + g[idx_g_R]    )  +
                                   (f[idx_f_L]     + f[idx_f_RR]   )*(g[idx_g_L]     + g[idx_g_RR]   )) +
                    half*(d_coef_c*(f[idx_f_LLL]   + f[idx_f_R]    )*(g[idx_g_LLL]   + g[idx_g_R]    )  +
                                   (f[idx_f_LL]    + f[idx_f_RR]   )*(g[idx_g_LL]    + g[idx_g_RR]   )  +
                                   (f[idx_f_L]     + f[idx_f_RRR]  )*(g[idx_g_L]     + g[idx_g_RRR]  )) +
                    half*(d_coef_d*(f[idx_f_LLLL]  + f[idx_f_R]    )*(g[idx_g_LLLL]  + g[idx_g_R]    )  +
                                   (f[idx_f_LLL]   + f[idx_f_RR]   )*(g[idx_g_LLL]   + g[idx_g_RR]   )  +
                                   (f[idx_f_LL]    + f[idx_f_RRR]  )*(g[idx_g_LL]    + g[idx_g_RRR]  )  +
                                   (f[idx_f_L]     + f[idx_f_RRRR] )*(g[idx_g_L]     + g[idx_g_RRRR] )) +
                    half*(d_coef_e*(f[idx_f_LLLLL] + f[idx_f_R]    )*(g[idx_g_LLLLL] + g[idx_g_R]    )  +
                                   (f[idx_f_LLLL]  + f[idx_f_RR]   )*(g[idx_g_LLLL]  + g[idx_g_RR]   )  +
                                   (f[idx_f_LLL]   + f[idx_f_RRR]  )*(g[idx_g_LLL]   + g[idx_g_RRR]  )  +
                                   (f[idx_f_LL]    + f[idx_f_RRRR] )*(g[idx_g_LL]    + g[idx_g_RRRR] )  +
                                   (f[idx_f_L]     + f[idx_f_RRRRR])*(g[idx_g_L]     + g[idx_g_RRRRR])));
            }
        }
        else if (d_stencil_width == 13)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LLLLLL = i - 6 + num_ghosts_0_f;
                const int idx_f_LLLLL  = i - 5 + num_ghosts_0_f;
                const int idx_f_LLLL   = i - 4 + num_ghosts_0_f;
                const int idx_f_LLL    = i - 3 + num_ghosts_0_f;
                const int idx_f_LL     = i - 2 + num_ghosts_0_f;
                const int idx_f_L      = i - 1 + num_ghosts_0_f;
                const int idx_f_R      = i     + num_ghosts_0_f;
                const int idx_f_RR     = i + 1 + num_ghosts_0_f;
                const int idx_f_RRR    = i + 2 + num_ghosts_0_f;
                const int idx_f_RRRR   = i + 3 + num_ghosts_0_f;
                const int idx_f_RRRRR  = i + 4 + num_ghosts_0_f;
                const int idx_f_RRRRRR = i + 5 + num_ghosts_0_f;
                
                const int idx_g_LLLLLL = i - 6 + num_ghosts_0_g;
                const int idx_g_LLLLL  = i - 5 + num_ghosts_0_g;
                const int idx_g_LLLL   = i - 4 + num_ghosts_0_g;
                const int idx_g_LLL    = i - 3 + num_ghosts_0_g;
                const int idx_g_LL     = i - 2 + num_ghosts_0_g;
                const int idx_g_L      = i - 1 + num_ghosts_0_g;
                const int idx_g_R      = i     + num_ghosts_0_g;
                const int idx_g_RR     = i + 1 + num_ghosts_0_g;
                const int idx_g_RRR    = i + 2 + num_ghosts_0_g;
                const int idx_g_RRRR   = i + 3 + num_ghosts_0_g;
                const int idx_g_RRRRR  = i + 4 + num_ghosts_0_g;
                const int idx_g_RRRRRR = i + 5 + num_ghosts_0_g;
                
                F_face_x[idx_face_x] += dt*(
                    half*(d_coef_a*(f[idx_f_L]      + f[idx_f_R]     )*(g[idx_g_L]      + g[idx_g_R]     )) +
                    half*(d_coef_b*(f[idx_f_LL]     + f[idx_f_R]     )*(g[idx_g_LL]     + g[idx_g_R]     )  +
                                   (f[idx_f_L]      + f[idx_f_RR]    )*(g[idx_g_L]      + g[idx_g_RR]    )) +
                    half*(d_coef_c*(f[idx_f_LLL]    + f[idx_f_R]     )*(g[idx_g_LLL]    + g[idx_g_R]     )  +
                                   (f[idx_f_LL]     + f[idx_f_RR]    )*(g[idx_g_LL]     + g[idx_g_RR]    )  +
                                   (f[idx_f_L]      + f[idx_f_RRR]   )*(g[idx_g_L]      + g[idx_g_RRR]   )) +
                    half*(d_coef_d*(f[idx_f_LLLL]   + f[idx_f_R]     )*(g[idx_g_LLLL]   + g[idx_g_R]     )  +
                                   (f[idx_f_LLL]    + f[idx_f_RR]    )*(g[idx_g_LLL]    + g[idx_g_RR]    )  +
                                   (f[idx_f_LL]     + f[idx_f_RRR]   )*(g[idx_g_LL]     + g[idx_g_RRR]   )  +
                                   (f[idx_f_L]      + f[idx_f_RRRR]  )*(g[idx_g_L]      + g[idx_g_RRRR]  )) +
                    half*(d_coef_e*(f[idx_f_LLLLL]  + f[idx_f_R]     )*(g[idx_g_LLLLL]  + g[idx_g_R]     )  +
                                   (f[idx_f_LLLL]   + f[idx_f_RR]    )*(g[idx_g_LLLL]   + g[idx_g_RR]    )  +
                                   (f[idx_f_LLL]    + f[idx_f_RRR]   )*(g[idx_g_LLL]    + g[idx_g_RRR]   )  +
                                   (f[idx_f_LL]     + f[idx_f_RRRR]  )*(g[idx_g_LL]     + g[idx_g_RRRR]  )  +
                                   (f[idx_f_L]      + f[idx_f_RRRRR] )*(g[idx_g_L]      + g[idx_g_RRRRR] )) +
                    half*(d_coef_f*(f[idx_f_LLLLLL] + f[idx_f_R]     )*(g[idx_g_LLLLLL] + g[idx_g_R]     )  +
                                   (f[idx_f_LLLLL]  + f[idx_f_RR]    )*(g[idx_g_LLLLL]  + g[idx_g_RR]    )  +
                                   (f[idx_f_LLLL]   + f[idx_f_RRR]   )*(g[idx_g_LLLL]   + g[idx_g_RRR]   )  +
                                   (f[idx_f_LLL]    + f[idx_f_RRRR]  )*(g[idx_g_LLL]    + g[idx_g_RRRR]  )  +
                                   (f[idx_f_LL]     + f[idx_f_RRRRR] )*(g[idx_g_LL]     + g[idx_g_RRRRR] )  +
                                   (f[idx_f_L]      + f[idx_f_RRRRRR])*(g[idx_g_L]      + g[idx_g_RRRRRR])));
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        /*
         * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
         */
        
        const int num_ghosts_0_f = num_ghosts_f[0];
        const int num_ghosts_1_f = num_ghosts_f[1];
        const int ghostcell_dim_0_f = ghostcell_dims_f[0];
        
        const int num_ghosts_0_g = num_ghosts_g[0];
        const int num_ghosts_1_g = num_ghosts_g[1];
        const int ghostcell_dim_0_g = ghostcell_dims_g[0];
        
        if (d_stencil_width == 3)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_L = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_L = (i - 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_R = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    F_face_x[idx_face_x] += dt*(
                        half*(d_coef_a*(f[idx_f_L] + f[idx_f_R])*(g[idx_g_L] + g[idx_g_R])));
                }
            }
        }
        else if (d_stencil_width == 5)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LL = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L  = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R  = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_LL = (i - 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_L  = (i - 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_R  = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RR = (i + 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    F_face_x[idx_face_x] += dt*(
                        half*(d_coef_a*(f[idx_f_L]  + f[idx_f_R] )*(g[idx_g_L]  + g[idx_g_R] )) +
                        half*(d_coef_b*(f[idx_f_LL] + f[idx_f_R] )*(g[idx_g_LL] + g[idx_g_R] )  +
                                       (f[idx_f_L]  + f[idx_f_RR])*(g[idx_g_L]  + g[idx_g_RR])));
                }
            }
        }
        else if (d_stencil_width == 7)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LLL = (i - 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LL  = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L   = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R   = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR  = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRR = (i + 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_LLL = (i - 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LL  = (i - 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_L   = (i - 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_R   = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RR  = (i + 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRR = (i + 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    F_face_x[idx_face_x] += dt*(
                        half*(d_coef_a*(f[idx_f_L]   + f[idx_f_R]  )*(g[idx_g_L]   + g[idx_g_R]  )) +
                        half*(d_coef_b*(f[idx_f_LL]  + f[idx_f_R]  )*(g[idx_g_LL]  + g[idx_g_R]  )  +
                                       (f[idx_f_L]   + f[idx_f_RR] )*(g[idx_g_L]   + g[idx_g_RR] )) +
                        half*(d_coef_c*(f[idx_f_LLL] + f[idx_f_R]  )*(g[idx_g_LLL] + g[idx_g_R]  )  +
                                       (f[idx_f_LL]  + f[idx_f_RR] )*(g[idx_g_LL]  + g[idx_g_RR] )  +
                                       (f[idx_f_L]   + f[idx_f_RRR])*(g[idx_g_L]   + g[idx_g_RRR])));
                }
            }
        }
        else if (d_stencil_width == 9)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LLLL = (i - 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLL  = (i - 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LL   = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L    = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R    = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR   = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRR  = (i + 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRR = (i + 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_LLLL = (i - 4 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LLL  = (i - 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LL   = (i - 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_L    = (i - 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_R    = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RR   = (i + 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRR  = (i + 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRRR = (i + 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    F_face_x[idx_face_x] += dt*(
                        half*(d_coef_a*(f[idx_f_L]    + f[idx_f_R]   )*(g[idx_g_L]    + g[idx_g_R]   )) +
                        half*(d_coef_b*(f[idx_f_LL]   + f[idx_f_R]   )*(g[idx_g_LL]   + g[idx_g_R]   )  +
                                       (f[idx_f_L]    + f[idx_f_RR]  )*(g[idx_g_L]    + g[idx_g_RR]  )) +
                        half*(d_coef_c*(f[idx_f_LLL]  + f[idx_f_R]   )*(g[idx_g_LLL]  + g[idx_g_R]   )  +
                                       (f[idx_f_LL]   + f[idx_f_RR]  )*(g[idx_g_LL]   + g[idx_g_RR]  )  +
                                       (f[idx_f_L]    + f[idx_f_RRR] )*(g[idx_g_L]    + g[idx_g_RRR] )) +
                        half*(d_coef_d*(f[idx_f_LLLL] + f[idx_f_R]   )*(g[idx_g_LLLL] + g[idx_g_R]   )  +
                                       (f[idx_f_LLL]  + f[idx_f_RR]  )*(g[idx_g_LLL]  + g[idx_g_RR]  )  +
                                       (f[idx_f_LL]   + f[idx_f_RRR] )*(g[idx_g_LL]   + g[idx_g_RRR] )  +
                                       (f[idx_f_L]    + f[idx_f_RRRR])*(g[idx_g_L]    + g[idx_g_RRRR])));
                }
            }
        }
        else if (d_stencil_width == 11)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LLLLL = (i - 5 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLLL  = (i - 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLL   = (i - 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LL    = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L     = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R     = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR    = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRR   = (i + 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRR  = (i + 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRRR = (i + 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_LLLLL = (i - 5 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LLLL  = (i - 4 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LLL   = (i - 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LL    = (i - 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_L     = (i - 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_R     = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RR    = (i + 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRR   = (i + 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRRR  = (i + 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRRRR = (i + 4 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    F_face_x[idx_face_x] += dt*(
                        half*(d_coef_a*(f[idx_f_L]     + f[idx_f_R]    )*(g[idx_g_L]     + g[idx_g_R]    )) +
                        half*(d_coef_b*(f[idx_f_LL]    + f[idx_f_R]    )*(g[idx_g_LL]    + g[idx_g_R]    )  +
                                       (f[idx_f_L]     + f[idx_f_RR]   )*(g[idx_g_L]     + g[idx_g_RR]   )) +
                        half*(d_coef_c*(f[idx_f_LLL]   + f[idx_f_R]    )*(g[idx_g_LLL]   + g[idx_g_R]    )  +
                                       (f[idx_f_LL]    + f[idx_f_RR]   )*(g[idx_g_LL]    + g[idx_g_RR]   )  +
                                       (f[idx_f_L]     + f[idx_f_RRR]  )*(g[idx_g_L]     + g[idx_g_RRR]  )) +
                        half*(d_coef_d*(f[idx_f_LLLL]  + f[idx_f_R]    )*(g[idx_g_LLLL]  + g[idx_g_R]    )  +
                                       (f[idx_f_LLL]   + f[idx_f_RR]   )*(g[idx_g_LLL]   + g[idx_g_RR]   )  +
                                       (f[idx_f_LL]    + f[idx_f_RRR]  )*(g[idx_g_LL]    + g[idx_g_RRR]  )  +
                                       (f[idx_f_L]     + f[idx_f_RRRR] )*(g[idx_g_L]     + g[idx_g_RRRR] )) +
                        half*(d_coef_e*(f[idx_f_LLLLL] + f[idx_f_R]    )*(g[idx_g_LLLLL] + g[idx_g_R]    )  +
                                       (f[idx_f_LLLL]  + f[idx_f_RR]   )*(g[idx_g_LLLL]  + g[idx_g_RR]   )  +
                                       (f[idx_f_LLL]   + f[idx_f_RRR]  )*(g[idx_g_LLL]   + g[idx_g_RRR]  )  +
                                       (f[idx_f_LL]    + f[idx_f_RRRR] )*(g[idx_g_LL]    + g[idx_g_RRRR] )  +
                                       (f[idx_f_L]     + f[idx_f_RRRRR])*(g[idx_g_L]     + g[idx_g_RRRRR])));
                }
            }
        }
        else if (d_stencil_width == 13)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LLLLLL = (i - 6 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLLLL  = (i - 5 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLLL   = (i - 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLL    = (i - 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LL     = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L      = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R      = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR     = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRR    = (i + 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRR   = (i + 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRRR  = (i + 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRRRR = (i + 5 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_LLLLLL = (i - 6 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LLLLL  = (i - 5 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LLLL   = (i - 4 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LLL    = (i - 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LL     = (i - 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_L      = (i - 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_R      = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RR     = (i + 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRR    = (i + 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRRR   = (i + 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRRRR  = (i + 4 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRRRRR = (i + 5 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    F_face_x[idx_face_x] += dt*(
                        half*(d_coef_a*(f[idx_f_L]      + f[idx_f_R]     )*(g[idx_g_L]      + g[idx_g_R]     )) +
                        half*(d_coef_b*(f[idx_f_LL]     + f[idx_f_R]     )*(g[idx_g_LL]     + g[idx_g_R]     )  +
                                       (f[idx_f_L]      + f[idx_f_RR]    )*(g[idx_g_L]      + g[idx_g_RR]    )) +
                        half*(d_coef_c*(f[idx_f_LLL]    + f[idx_f_R]     )*(g[idx_g_LLL]    + g[idx_g_R]     )  +
                                       (f[idx_f_LL]     + f[idx_f_RR]    )*(g[idx_g_LL]     + g[idx_g_RR]    )  +
                                       (f[idx_f_L]      + f[idx_f_RRR]   )*(g[idx_g_L]      + g[idx_g_RRR]   )) +
                        half*(d_coef_d*(f[idx_f_LLLL]   + f[idx_f_R]     )*(g[idx_g_LLLL]   + g[idx_g_R]     )  +
                                       (f[idx_f_LLL]    + f[idx_f_RR]    )*(g[idx_g_LLL]    + g[idx_g_RR]    )  +
                                       (f[idx_f_LL]     + f[idx_f_RRR]   )*(g[idx_g_LL]     + g[idx_g_RRR]   )  +
                                       (f[idx_f_L]      + f[idx_f_RRRR]  )*(g[idx_g_L]      + g[idx_g_RRRR]  )) +
                        half*(d_coef_e*(f[idx_f_LLLLL]  + f[idx_f_R]     )*(g[idx_g_LLLLL]  + g[idx_g_R]     )  +
                                       (f[idx_f_LLLL]   + f[idx_f_RR]    )*(g[idx_g_LLLL]   + g[idx_g_RR]    )  +
                                       (f[idx_f_LLL]    + f[idx_f_RRR]   )*(g[idx_g_LLL]    + g[idx_g_RRR]   )  +
                                       (f[idx_f_LL]     + f[idx_f_RRRR]  )*(g[idx_g_LL]     + g[idx_g_RRRR]  )  +
                                       (f[idx_f_L]      + f[idx_f_RRRRR] )*(g[idx_g_L]      + g[idx_g_RRRRR] )) +
                        half*(d_coef_f*(f[idx_f_LLLLLL] + f[idx_f_R]     )*(g[idx_g_LLLLLL] + g[idx_g_R]     )  +
                                       (f[idx_f_LLLLL]  + f[idx_f_RR]    )*(g[idx_g_LLLLL]  + g[idx_g_RR]    )  +
                                       (f[idx_f_LLLL]   + f[idx_f_RRR]   )*(g[idx_g_LLLL]   + g[idx_g_RRR]   )  +
                                       (f[idx_f_LLL]    + f[idx_f_RRRR]  )*(g[idx_g_LLL]    + g[idx_g_RRRR]  )  +
                                       (f[idx_f_LL]     + f[idx_f_RRRRR] )*(g[idx_g_LL]     + g[idx_g_RRRRR] )  +
                                       (f[idx_f_L]      + f[idx_f_RRRRRR])*(g[idx_g_L]      + g[idx_g_RRRRRR])));
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        /*
         * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
         */
        
        const int num_ghosts_0_f = num_ghosts_f[0];
        const int num_ghosts_1_f = num_ghosts_f[1];
        const int num_ghosts_2_f = num_ghosts_f[2];
        const int ghostcell_dim_0_f = ghostcell_dims_f[0];
        const int ghostcell_dim_1_f = ghostcell_dims_f[1];
        
        const int num_ghosts_0_g = num_ghosts_g[0];
        const int num_ghosts_1_g = num_ghosts_g[1];
        const int num_ghosts_2_g = num_ghosts_g[2];
        const int ghostcell_dim_0_g = ghostcell_dims_g[0];
        const int ghostcell_dim_1_g = ghostcell_dims_g[1];
        
        if (d_stencil_width == 3)
        {
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
                        
                        const int idx_f_L = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_L = (i - 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_R = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_x[idx_face_x] += dt*(
                            half*(d_coef_a*(f[idx_f_L] + f[idx_f_R])*(g[idx_g_L] + g[idx_g_R])));
                    }
                }
            }
        }
        else if (d_stencil_width == 5)
        {
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
                        
                        const int idx_f_LL = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L  = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_LL = (i - 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_L  = (i - 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_R  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RR = (i + 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_x[idx_face_x] += dt*(
                            half*(d_coef_a*(f[idx_f_L]  + f[idx_f_R] )*(g[idx_g_L]  + g[idx_g_R] )) +
                            half*(d_coef_b*(f[idx_f_LL] + f[idx_f_R] )*(g[idx_g_LL] + g[idx_g_R] )  +
                                           (f[idx_f_L]  + f[idx_f_RR])*(g[idx_g_L]  + g[idx_g_RR])));
                    }
                }
            }
        }
        else if (d_stencil_width == 7)
        {
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
                        
                        const int idx_f_LLL = (i - 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LL  = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L   = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR  = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRR = (i + 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_LLL = (i - 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LL  = (i - 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_L   = (i - 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_R   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RR  = (i + 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRR = (i + 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_x[idx_face_x] += dt*(
                            half*(d_coef_a*(f[idx_f_L]   + f[idx_f_R]  )*(g[idx_g_L]   + g[idx_g_R]  )) +
                            half*(d_coef_b*(f[idx_f_LL]  + f[idx_f_R]  )*(g[idx_g_LL]  + g[idx_g_R]  )  +
                                           (f[idx_f_L]   + f[idx_f_RR] )*(g[idx_g_L]   + g[idx_g_RR] )) +
                            half*(d_coef_c*(f[idx_f_LLL] + f[idx_f_R]  )*(g[idx_g_LLL] + g[idx_g_R]  )  +
                                           (f[idx_f_LL]  + f[idx_f_RR] )*(g[idx_g_LL]  + g[idx_g_RR] )  +
                                           (f[idx_f_L]   + f[idx_f_RRR])*(g[idx_g_L]   + g[idx_g_RRR])));
                    }
                }
            }
        }
        else if (d_stencil_width == 9)
        {
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
                        
                        const int idx_f_LLLL = (i - 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLL  = (i - 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LL   = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L    = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR   = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRR  = (i + 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRR = (i + 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_LLLL = (i - 4 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LLL  = (i - 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LL   = (i - 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_L    = (i - 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_R    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RR   = (i + 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRR  = (i + 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRRR = (i + 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_x[idx_face_x] += dt*(
                            half*(d_coef_a*(f[idx_f_L]    + f[idx_f_R]   )*(g[idx_g_L]    + g[idx_g_R]   )) +
                            half*(d_coef_b*(f[idx_f_LL]   + f[idx_f_R]   )*(g[idx_g_LL]   + g[idx_g_R]   )  +
                                           (f[idx_f_L]    + f[idx_f_RR]  )*(g[idx_g_L]    + g[idx_g_RR]  )) +
                            half*(d_coef_c*(f[idx_f_LLL]  + f[idx_f_R]   )*(g[idx_g_LLL]  + g[idx_g_R]   )  +
                                           (f[idx_f_LL]   + f[idx_f_RR]  )*(g[idx_g_LL]   + g[idx_g_RR]  )  +
                                           (f[idx_f_L]    + f[idx_f_RRR] )*(g[idx_g_L]    + g[idx_g_RRR] )) +
                            half*(d_coef_d*(f[idx_f_LLLL] + f[idx_f_R]   )*(g[idx_g_LLLL] + g[idx_g_R]   )  +
                                           (f[idx_f_LLL]  + f[idx_f_RR]  )*(g[idx_g_LLL]  + g[idx_g_RR]  )  +
                                           (f[idx_f_LL]   + f[idx_f_RRR] )*(g[idx_g_LL]   + g[idx_g_RRR] )  +
                                           (f[idx_f_L]    + f[idx_f_RRRR])*(g[idx_g_L]    + g[idx_g_RRRR])));
                    }
                }
            }
        }
        else if (d_stencil_width == 11)
        {
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
                        
                        const int idx_f_LLLLL = (i - 5 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLLL  = (i - 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLL   = (i - 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LL    = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L     = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR    = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRR   = (i + 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRR  = (i + 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRRR = (i + 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_LLLLL = (i - 5 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LLLL  = (i - 4 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LLL   = (i - 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LL    = (i - 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_L     = (i - 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_R     = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RR    = (i + 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRR   = (i + 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRRR  = (i + 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRRRR = (i + 4 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_x[idx_face_x] += dt*(
                            half*(d_coef_a*(f[idx_f_L]     + f[idx_f_R]    )*(g[idx_g_L]     + g[idx_g_R]    )) +
                            half*(d_coef_b*(f[idx_f_LL]    + f[idx_f_R]    )*(g[idx_g_LL]    + g[idx_g_R]    )  +
                                           (f[idx_f_L]     + f[idx_f_RR]   )*(g[idx_g_L]     + g[idx_g_RR]   )) +
                            half*(d_coef_c*(f[idx_f_LLL]   + f[idx_f_R]    )*(g[idx_g_LLL]   + g[idx_g_R]    )  +
                                           (f[idx_f_LL]    + f[idx_f_RR]   )*(g[idx_g_LL]    + g[idx_g_RR]   )  +
                                           (f[idx_f_L]     + f[idx_f_RRR]  )*(g[idx_g_L]     + g[idx_g_RRR]  )) +
                            half*(d_coef_d*(f[idx_f_LLLL]  + f[idx_f_R]    )*(g[idx_g_LLLL]  + g[idx_g_R]    )  +
                                           (f[idx_f_LLL]   + f[idx_f_RR]   )*(g[idx_g_LLL]   + g[idx_g_RR]   )  +
                                           (f[idx_f_LL]    + f[idx_f_RRR]  )*(g[idx_g_LL]    + g[idx_g_RRR]  )  +
                                           (f[idx_f_L]     + f[idx_f_RRRR] )*(g[idx_g_L]     + g[idx_g_RRRR] )) +
                            half*(d_coef_e*(f[idx_f_LLLLL] + f[idx_f_R]    )*(g[idx_g_LLLLL] + g[idx_g_R]    )  +
                                           (f[idx_f_LLLL]  + f[idx_f_RR]   )*(g[idx_g_LLLL]  + g[idx_g_RR]   )  +
                                           (f[idx_f_LLL]   + f[idx_f_RRR]  )*(g[idx_g_LLL]   + g[idx_g_RRR]  )  +
                                           (f[idx_f_LL]    + f[idx_f_RRRR] )*(g[idx_g_LL]    + g[idx_g_RRRR] )  +
                                           (f[idx_f_L]     + f[idx_f_RRRRR])*(g[idx_g_L]     + g[idx_g_RRRRR])));
                    }
                }
            }
        }
        else if (d_stencil_width == 13)
        {
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
                        
                        const int idx_f_LLLLLL = (i - 6 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLLLL  = (i - 5 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLLL   = (i - 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLL    = (i - 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LL     = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L      = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R      = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR     = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRR    = (i + 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRR   = (i + 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRRR  = (i + 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRRRR = (i + 5 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_LLLLLL = (i - 6 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LLLLL  = (i - 5 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LLLL   = (i - 4 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LLL    = (i - 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LL     = (i - 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_L      = (i - 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_R      = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RR     = (i + 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRR    = (i + 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRRR   = (i + 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRRRR  = (i + 4 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRRRRR = (i + 5 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_x[idx_face_x] += dt*(
                            half*(d_coef_a*(f[idx_f_L]      + f[idx_f_R]     )*(g[idx_g_L]      + g[idx_g_R]     )) +
                            half*(d_coef_b*(f[idx_f_LL]     + f[idx_f_R]     )*(g[idx_g_LL]     + g[idx_g_R]     )  +
                                           (f[idx_f_L]      + f[idx_f_RR]    )*(g[idx_g_L]      + g[idx_g_RR]    )) +
                            half*(d_coef_c*(f[idx_f_LLL]    + f[idx_f_R]     )*(g[idx_g_LLL]    + g[idx_g_R]     )  +
                                           (f[idx_f_LL]     + f[idx_f_RR]    )*(g[idx_g_LL]     + g[idx_g_RR]    )  +
                                           (f[idx_f_L]      + f[idx_f_RRR]   )*(g[idx_g_L]      + g[idx_g_RRR]   )) +
                            half*(d_coef_d*(f[idx_f_LLLL]   + f[idx_f_R]     )*(g[idx_g_LLLL]   + g[idx_g_R]     )  +
                                           (f[idx_f_LLL]    + f[idx_f_RR]    )*(g[idx_g_LLL]    + g[idx_g_RR]    )  +
                                           (f[idx_f_LL]     + f[idx_f_RRR]   )*(g[idx_g_LL]     + g[idx_g_RRR]   )  +
                                           (f[idx_f_L]      + f[idx_f_RRRR]  )*(g[idx_g_L]      + g[idx_g_RRRR]  )) +
                            half*(d_coef_e*(f[idx_f_LLLLL]  + f[idx_f_R]     )*(g[idx_g_LLLLL]  + g[idx_g_R]     )  +
                                           (f[idx_f_LLLL]   + f[idx_f_RR]    )*(g[idx_g_LLLL]   + g[idx_g_RR]    )  +
                                           (f[idx_f_LLL]    + f[idx_f_RRR]   )*(g[idx_g_LLL]    + g[idx_g_RRR]   )  +
                                           (f[idx_f_LL]     + f[idx_f_RRRR]  )*(g[idx_g_LL]     + g[idx_g_RRRR]  )  +
                                           (f[idx_f_L]      + f[idx_f_RRRRR] )*(g[idx_g_L]      + g[idx_g_RRRRR] )) +
                            half*(d_coef_f*(f[idx_f_LLLLLL] + f[idx_f_R]     )*(g[idx_g_LLLLLL] + g[idx_g_R]     )  +
                                           (f[idx_f_LLLLL]  + f[idx_f_RR]    )*(g[idx_g_LLLLL]  + g[idx_g_RR]    )  +
                                           (f[idx_f_LLLL]   + f[idx_f_RRR]   )*(g[idx_g_LLLL]   + g[idx_g_RRR]   )  +
                                           (f[idx_f_LLL]    + f[idx_f_RRRR]  )*(g[idx_g_LLL]    + g[idx_g_RRRR]  )  +
                                           (f[idx_f_LL]     + f[idx_f_RRRRR] )*(g[idx_g_LL]     + g[idx_g_RRRRR] )  +
                                           (f[idx_f_L]      + f[idx_f_RRRRRR])*(g[idx_g_L]      + g[idx_g_RRRRRR])));
                    }
                }
            }
        }
    }
}


/*
 * Add cubic term to convective flux in x-direction.
 */
void
ConvectiveFluxReconstructorKEP::addCubicTermToConvectiveFluxX(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_f,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_g,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_h,
    const int component_idx_flux,
    const int component_idx_f,
    const int component_idx_g,
    const int component_idx_h,
    const double dt) const
{
    const hier::Box interior_box = data_convective_flux->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(data_f->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_g->getBox().isSpatiallyEqual(interior_box));
        TBOX_ASSERT(data_h->getBox().isSpatiallyEqual(interior_box));
#endif
    
    // Get the dimensions of the ghost cell boxes.
    const hier::Box ghost_box_f = data_f->getGhostBox();
    const hier::IntVector ghostcell_dims_f = ghost_box_f.numberCells();
    const hier::IntVector num_ghosts_f = data_f->getGhostCellWidth();
    
    const hier::Box ghost_box_g = data_g->getGhostBox();
    const hier::IntVector ghostcell_dims_g = ghost_box_g.numberCells();
    const hier::IntVector num_ghosts_g = data_g->getGhostCellWidth();
    
    const hier::Box ghost_box_h = data_h->getGhostBox();
    const hier::IntVector ghostcell_dims_h = ghost_box_h.numberCells();
    const hier::IntVector num_ghosts_h = data_h->getGhostCellWidth();
    
    /*
     * Get the pointers to the data.
     */
    
    double* F_face_x = data_convective_flux->getPointer(0, component_idx_flux);
    double* f        = data_f->getPointer(component_idx_f);
    double* g        = data_g->getPointer(component_idx_g);
    double* h        = data_h->getPointer(component_idx_h);
    
    const double quarter = double(1)/double(4);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Get the numbers of ghost cells.
         */
        const int num_ghosts_0_f = num_ghosts_f[0];
        const int num_ghosts_0_g = num_ghosts_g[0];
        const int num_ghosts_0_h = num_ghosts_h[0];
        
        if (d_stencil_width == 3)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_L = i - 1 + num_ghosts_0_f;
                const int idx_f_R = i     + num_ghosts_0_f;
                
                const int idx_g_L = i - 1 + num_ghosts_0_g;
                const int idx_g_R = i     + num_ghosts_0_g;
                
                const int idx_h_L = i - 1 + num_ghosts_0_h;
                const int idx_h_R = i     + num_ghosts_0_h;
                
                F_face_x[idx_face_x] += dt*(
                    quarter*(d_coef_a*(f[idx_f_L] + f[idx_f_R])*(g[idx_g_L] + g[idx_g_R])*(h[idx_h_L] + h[idx_h_R])));
            }
        }
        else if (d_stencil_width == 5)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LL = i - 2 + num_ghosts_0_f;
                const int idx_f_L  = i - 1 + num_ghosts_0_f;
                const int idx_f_R  = i     + num_ghosts_0_f;
                const int idx_f_RR = i + 1 + num_ghosts_0_f;
                
                const int idx_g_LL = i - 2 + num_ghosts_0_g;
                const int idx_g_L  = i - 1 + num_ghosts_0_g;
                const int idx_g_R  = i     + num_ghosts_0_g;
                const int idx_g_RR = i + 1 + num_ghosts_0_g;
                
                const int idx_h_LL = i - 2 + num_ghosts_0_h;
                const int idx_h_L  = i - 1 + num_ghosts_0_h;
                const int idx_h_R  = i     + num_ghosts_0_h;
                const int idx_h_RR = i + 1 + num_ghosts_0_h;
                
                F_face_x[idx_face_x] += dt*(
                    quarter*(d_coef_a*(f[idx_f_L]  + f[idx_f_R] )*(g[idx_g_L]  + g[idx_g_R] )*(h[idx_h_L]  + h[idx_h_R] )) +
                    quarter*(d_coef_b*(f[idx_f_LL] + f[idx_f_R] )*(g[idx_g_LL] + g[idx_g_R] )*(h[idx_h_LL] + h[idx_h_R] )  +
                                      (f[idx_f_L]  + f[idx_f_RR])*(g[idx_g_L]  + g[idx_g_RR])*(h[idx_h_L]  + h[idx_h_RR])));
            }
        }
        else if (d_stencil_width == 7)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LLL = i - 3 + num_ghosts_0_f;
                const int idx_f_LL  = i - 2 + num_ghosts_0_f;
                const int idx_f_L   = i - 1 + num_ghosts_0_f;
                const int idx_f_R   = i     + num_ghosts_0_f;
                const int idx_f_RR  = i + 1 + num_ghosts_0_f;
                const int idx_f_RRR = i + 2 + num_ghosts_0_f;
                
                const int idx_g_LLL = i - 3 + num_ghosts_0_g;
                const int idx_g_LL  = i - 2 + num_ghosts_0_g;
                const int idx_g_L   = i - 1 + num_ghosts_0_g;
                const int idx_g_R   = i     + num_ghosts_0_g;
                const int idx_g_RR  = i + 1 + num_ghosts_0_g;
                const int idx_g_RRR = i + 2 + num_ghosts_0_g;
                
                const int idx_h_LLL = i - 3 + num_ghosts_0_h;
                const int idx_h_LL  = i - 2 + num_ghosts_0_h;
                const int idx_h_L   = i - 1 + num_ghosts_0_h;
                const int idx_h_R   = i     + num_ghosts_0_h;
                const int idx_h_RR  = i + 1 + num_ghosts_0_h;
                const int idx_h_RRR = i + 2 + num_ghosts_0_h;
                
                F_face_x[idx_face_x] += dt*(
                    quarter*(d_coef_a*(f[idx_f_L]   + f[idx_f_R]  )*(g[idx_g_L]   + g[idx_g_R]  )*(h[idx_h_L]   + h[idx_h_R]  )) +
                    quarter*(d_coef_b*(f[idx_f_LL]  + f[idx_f_R]  )*(g[idx_g_LL]  + g[idx_g_R]  )*(h[idx_h_LL]  + h[idx_h_R]  )  +
                                      (f[idx_f_L]   + f[idx_f_RR] )*(g[idx_g_L]   + g[idx_g_RR] )*(h[idx_h_L]   + h[idx_h_RR] )) +
                    quarter*(d_coef_c*(f[idx_f_LLL] + f[idx_f_R]  )*(g[idx_g_LLL] + g[idx_g_R]  )*(h[idx_h_LLL] + h[idx_h_R]  )  +
                                      (f[idx_f_LL]  + f[idx_f_RR] )*(g[idx_g_LL]  + g[idx_g_RR] )*(h[idx_h_LL]  + h[idx_h_RR] )  +
                                      (f[idx_f_L]   + f[idx_f_RRR])*(g[idx_g_L]   + g[idx_g_RRR])*(h[idx_h_L]   + h[idx_h_RRR])));
            }
        }
        else if (d_stencil_width == 9)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LLLL = i - 4 + num_ghosts_0_f;
                const int idx_f_LLL  = i - 3 + num_ghosts_0_f;
                const int idx_f_LL   = i - 2 + num_ghosts_0_f;
                const int idx_f_L    = i - 1 + num_ghosts_0_f;
                const int idx_f_R    = i     + num_ghosts_0_f;
                const int idx_f_RR   = i + 1 + num_ghosts_0_f;
                const int idx_f_RRR  = i + 2 + num_ghosts_0_f;
                const int idx_f_RRRR = i + 3 + num_ghosts_0_f;
                
                const int idx_g_LLLL = i - 4 + num_ghosts_0_g;
                const int idx_g_LLL  = i - 3 + num_ghosts_0_g;
                const int idx_g_LL   = i - 2 + num_ghosts_0_g;
                const int idx_g_L    = i - 1 + num_ghosts_0_g;
                const int idx_g_R    = i     + num_ghosts_0_g;
                const int idx_g_RR   = i + 1 + num_ghosts_0_g;
                const int idx_g_RRR  = i + 2 + num_ghosts_0_g;
                const int idx_g_RRRR = i + 3 + num_ghosts_0_g;
                
                const int idx_h_LLLL = i - 4 + num_ghosts_0_h;
                const int idx_h_LLL  = i - 3 + num_ghosts_0_h;
                const int idx_h_LL   = i - 2 + num_ghosts_0_h;
                const int idx_h_L    = i - 1 + num_ghosts_0_h;
                const int idx_h_R    = i     + num_ghosts_0_h;
                const int idx_h_RR   = i + 1 + num_ghosts_0_h;
                const int idx_h_RRR  = i + 2 + num_ghosts_0_h;
                const int idx_h_RRRR = i + 3 + num_ghosts_0_h;
                
                F_face_x[idx_face_x] += dt*(
                    quarter*(d_coef_a*(f[idx_f_L]    + f[idx_f_R]   )*(g[idx_g_L]    + g[idx_g_R]   )*(h[idx_h_L]    + h[idx_h_R]   )) +
                    quarter*(d_coef_b*(f[idx_f_LL]   + f[idx_f_R]   )*(g[idx_g_LL]   + g[idx_g_R]   )*(h[idx_h_LL]   + h[idx_h_R]   )  +
                                      (f[idx_f_L]    + f[idx_f_RR]  )*(g[idx_g_L]    + g[idx_g_RR]  )*(h[idx_h_L]    + h[idx_h_RR]  )) +
                    quarter*(d_coef_c*(f[idx_f_LLL]  + f[idx_f_R]   )*(g[idx_g_LLL]  + g[idx_g_R]   )*(h[idx_h_LLL]  + h[idx_h_R]   )  +
                                      (f[idx_f_LL]   + f[idx_f_RR]  )*(g[idx_g_LL]   + g[idx_g_RR]  )*(h[idx_h_LL]   + h[idx_h_RR]  )  +
                                      (f[idx_f_L]    + f[idx_f_RRR] )*(g[idx_g_L]    + g[idx_g_RRR] )*(h[idx_h_L]    + h[idx_h_RRR] )) +
                    quarter*(d_coef_d*(f[idx_f_LLLL] + f[idx_f_R]   )*(g[idx_g_LLLL] + g[idx_g_R]   )*(h[idx_h_LLLL] + h[idx_h_R]   )  +
                                      (f[idx_f_LLL]  + f[idx_f_RR]  )*(g[idx_g_LLL]  + g[idx_g_RR]  )*(h[idx_h_LLL]  + h[idx_h_RR]  )  +
                                      (f[idx_f_LL]   + f[idx_f_RRR] )*(g[idx_g_LL]   + g[idx_g_RRR] )*(h[idx_h_LL]   + h[idx_h_RRR] )  +
                                      (f[idx_f_L]    + f[idx_f_RRRR])*(g[idx_g_L]    + g[idx_g_RRRR])*(h[idx_h_L]    + h[idx_h_RRRR])));
            }
        }
        else if (d_stencil_width == 11)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LLLLL = i - 5 + num_ghosts_0_f;
                const int idx_f_LLLL  = i - 4 + num_ghosts_0_f;
                const int idx_f_LLL   = i - 3 + num_ghosts_0_f;
                const int idx_f_LL    = i - 2 + num_ghosts_0_f;
                const int idx_f_L     = i - 1 + num_ghosts_0_f;
                const int idx_f_R     = i     + num_ghosts_0_f;
                const int idx_f_RR    = i + 1 + num_ghosts_0_f;
                const int idx_f_RRR   = i + 2 + num_ghosts_0_f;
                const int idx_f_RRRR  = i + 3 + num_ghosts_0_f;
                const int idx_f_RRRRR = i + 4 + num_ghosts_0_f;
                
                const int idx_g_LLLLL = i - 5 + num_ghosts_0_g;
                const int idx_g_LLLL  = i - 4 + num_ghosts_0_g;
                const int idx_g_LLL   = i - 3 + num_ghosts_0_g;
                const int idx_g_LL    = i - 2 + num_ghosts_0_g;
                const int idx_g_L     = i - 1 + num_ghosts_0_g;
                const int idx_g_R     = i     + num_ghosts_0_g;
                const int idx_g_RR    = i + 1 + num_ghosts_0_g;
                const int idx_g_RRR   = i + 2 + num_ghosts_0_g;
                const int idx_g_RRRR  = i + 3 + num_ghosts_0_g;
                const int idx_g_RRRRR = i + 4 + num_ghosts_0_g;
                
                const int idx_h_LLLLL = i - 5 + num_ghosts_0_h;
                const int idx_h_LLLL  = i - 4 + num_ghosts_0_h;
                const int idx_h_LLL   = i - 3 + num_ghosts_0_h;
                const int idx_h_LL    = i - 2 + num_ghosts_0_h;
                const int idx_h_L     = i - 1 + num_ghosts_0_h;
                const int idx_h_R     = i     + num_ghosts_0_h;
                const int idx_h_RR    = i + 1 + num_ghosts_0_h;
                const int idx_h_RRR   = i + 2 + num_ghosts_0_h;
                const int idx_h_RRRR  = i + 3 + num_ghosts_0_h;
                const int idx_h_RRRRR = i + 4 + num_ghosts_0_h;
                
                F_face_x[idx_face_x] += dt*(
                    quarter*(d_coef_a*(f[idx_f_L]     + f[idx_f_R]    )*(g[idx_g_L]     + g[idx_g_R]    )*(h[idx_h_L]     + h[idx_h_R]    )) +
                    quarter*(d_coef_b*(f[idx_f_LL]    + f[idx_f_R]    )*(g[idx_g_LL]    + g[idx_g_R]    )*(h[idx_h_LL]    + h[idx_h_R]    )  +
                                      (f[idx_f_L]     + f[idx_f_RR]   )*(g[idx_g_L]     + g[idx_g_RR]   )*(h[idx_h_L]     + h[idx_h_RR]   )) +
                    quarter*(d_coef_c*(f[idx_f_LLL]   + f[idx_f_R]    )*(g[idx_g_LLL]   + g[idx_g_R]    )*(h[idx_h_LLL]   + h[idx_h_R]    )  +
                                      (f[idx_f_LL]    + f[idx_f_RR]   )*(g[idx_g_LL]    + g[idx_g_RR]   )*(h[idx_h_LL]    + h[idx_h_RR]   )  +
                                      (f[idx_f_L]     + f[idx_f_RRR]  )*(g[idx_g_L]     + g[idx_g_RRR]  )*(h[idx_h_L]     + h[idx_h_RRR]  )) +
                    quarter*(d_coef_d*(f[idx_f_LLLL]  + f[idx_f_R]    )*(g[idx_g_LLLL]  + g[idx_g_R]    )*(h[idx_h_LLLL]  + h[idx_h_R]    )  +
                                      (f[idx_f_LLL]   + f[idx_f_RR]   )*(g[idx_g_LLL]   + g[idx_g_RR]   )*(h[idx_h_LLL]   + h[idx_h_RR]   )  +
                                      (f[idx_f_LL]    + f[idx_f_RRR]  )*(g[idx_g_LL]    + g[idx_g_RRR]  )*(h[idx_h_LL]    + h[idx_h_RRR]  )  +
                                      (f[idx_f_L]     + f[idx_f_RRRR] )*(g[idx_g_L]     + g[idx_g_RRRR] )*(h[idx_h_L]     + h[idx_h_RRRR] )) +
                    quarter*(d_coef_e*(f[idx_f_LLLLL] + f[idx_f_R]    )*(g[idx_g_LLLLL] + g[idx_g_R]    )*(h[idx_h_LLLLL] + h[idx_h_R]    )  +
                                      (f[idx_f_LLLL]  + f[idx_f_RR]   )*(g[idx_g_LLLL]  + g[idx_g_RR]   )*(h[idx_h_LLLL]  + h[idx_h_RR]   )  +
                                      (f[idx_f_LLL]   + f[idx_f_RRR]  )*(g[idx_g_LLL]   + g[idx_g_RRR]  )*(h[idx_h_LLL]   + h[idx_h_RRR]  )  +
                                      (f[idx_f_LL]    + f[idx_f_RRRR] )*(g[idx_g_LL]    + g[idx_g_RRRR] )*(h[idx_h_LL]    + h[idx_h_RRRR] )  +
                                      (f[idx_f_L]     + f[idx_f_RRRRR])*(g[idx_g_L]     + g[idx_g_RRRRR])*(h[idx_h_L]     + h[idx_h_RRRRR])));
            }
        }
        else if (d_stencil_width == 13)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = 0; i < interior_dim_0 + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                
                const int idx_f_LLLLLL = i - 6 + num_ghosts_0_f;
                const int idx_f_LLLLL  = i - 5 + num_ghosts_0_f;
                const int idx_f_LLLL   = i - 4 + num_ghosts_0_f;
                const int idx_f_LLL    = i - 3 + num_ghosts_0_f;
                const int idx_f_LL     = i - 2 + num_ghosts_0_f;
                const int idx_f_L      = i - 1 + num_ghosts_0_f;
                const int idx_f_R      = i     + num_ghosts_0_f;
                const int idx_f_RR     = i + 1 + num_ghosts_0_f;
                const int idx_f_RRR    = i + 2 + num_ghosts_0_f;
                const int idx_f_RRRR   = i + 3 + num_ghosts_0_f;
                const int idx_f_RRRRR  = i + 4 + num_ghosts_0_f;
                const int idx_f_RRRRRR = i + 5 + num_ghosts_0_f;
                
                const int idx_g_LLLLLL = i - 6 + num_ghosts_0_g;
                const int idx_g_LLLLL  = i - 5 + num_ghosts_0_g;
                const int idx_g_LLLL   = i - 4 + num_ghosts_0_g;
                const int idx_g_LLL    = i - 3 + num_ghosts_0_g;
                const int idx_g_LL     = i - 2 + num_ghosts_0_g;
                const int idx_g_L      = i - 1 + num_ghosts_0_g;
                const int idx_g_R      = i     + num_ghosts_0_g;
                const int idx_g_RR     = i + 1 + num_ghosts_0_g;
                const int idx_g_RRR    = i + 2 + num_ghosts_0_g;
                const int idx_g_RRRR   = i + 3 + num_ghosts_0_g;
                const int idx_g_RRRRR  = i + 4 + num_ghosts_0_g;
                const int idx_g_RRRRRR = i + 5 + num_ghosts_0_g;
                
                const int idx_h_LLLLLL = i - 6 + num_ghosts_0_h;
                const int idx_h_LLLLL  = i - 5 + num_ghosts_0_h;
                const int idx_h_LLLL   = i - 4 + num_ghosts_0_h;
                const int idx_h_LLL    = i - 3 + num_ghosts_0_h;
                const int idx_h_LL     = i - 2 + num_ghosts_0_h;
                const int idx_h_L      = i - 1 + num_ghosts_0_h;
                const int idx_h_R      = i     + num_ghosts_0_h;
                const int idx_h_RR     = i + 1 + num_ghosts_0_h;
                const int idx_h_RRR    = i + 2 + num_ghosts_0_h;
                const int idx_h_RRRR   = i + 3 + num_ghosts_0_h;
                const int idx_h_RRRRR  = i + 4 + num_ghosts_0_h;
                const int idx_h_RRRRRR = i + 5 + num_ghosts_0_h;
                
                F_face_x[idx_face_x] += dt*(
                    quarter*(d_coef_a*(f[idx_f_L]      + f[idx_f_R]     )*(g[idx_g_L]      + g[idx_g_R]     )*(h[idx_h_L]      + h[idx_h_R]     )) +
                    quarter*(d_coef_b*(f[idx_f_LL]     + f[idx_f_R]     )*(g[idx_g_LL]     + g[idx_g_R]     )*(h[idx_h_LL]     + h[idx_h_R]     )  +
                                      (f[idx_f_L]      + f[idx_f_RR]    )*(g[idx_g_L]      + g[idx_g_RR]    )*(h[idx_h_L]      + h[idx_h_RR]    )) +
                    quarter*(d_coef_c*(f[idx_f_LLL]    + f[idx_f_R]     )*(g[idx_g_LLL]    + g[idx_g_R]     )*(h[idx_h_LLL]    + h[idx_h_R]     )  +
                                      (f[idx_f_LL]     + f[idx_f_RR]    )*(g[idx_g_LL]     + g[idx_g_RR]    )*(h[idx_h_LL]     + h[idx_h_RR]    )  +
                                      (f[idx_f_L]      + f[idx_f_RRR]   )*(g[idx_g_L]      + g[idx_g_RRR]   )*(h[idx_h_L]      + h[idx_h_RRR]   )) +
                    quarter*(d_coef_d*(f[idx_f_LLLL]   + f[idx_f_R]     )*(g[idx_g_LLLL]   + g[idx_g_R]     )*(h[idx_h_LLLL]   + h[idx_h_R]     )  +
                                      (f[idx_f_LLL]    + f[idx_f_RR]    )*(g[idx_g_LLL]    + g[idx_g_RR]    )*(h[idx_h_LLL]    + h[idx_h_RR]    )  +
                                      (f[idx_f_LL]     + f[idx_f_RRR]   )*(g[idx_g_LL]     + g[idx_g_RRR]   )*(h[idx_h_LL]     + h[idx_h_RRR]   )  +
                                      (f[idx_f_L]      + f[idx_f_RRRR]  )*(g[idx_g_L]      + g[idx_g_RRRR]  )*(h[idx_h_L]      + h[idx_h_RRRR]  )) +
                    quarter*(d_coef_e*(f[idx_f_LLLLL]  + f[idx_f_R]     )*(g[idx_g_LLLLL]  + g[idx_g_R]     )*(h[idx_h_LLLLL]  + h[idx_h_R]     )  +
                                      (f[idx_f_LLLL]   + f[idx_f_RR]    )*(g[idx_g_LLLL]   + g[idx_g_RR]    )*(h[idx_h_LLLL]   + h[idx_h_RR]    )  +
                                      (f[idx_f_LLL]    + f[idx_f_RRR]   )*(g[idx_g_LLL]    + g[idx_g_RRR]   )*(h[idx_h_LLL]    + h[idx_h_RRR]   )  +
                                      (f[idx_f_LL]     + f[idx_f_RRRR]  )*(g[idx_g_LL]     + g[idx_g_RRRR]  )*(h[idx_h_LL]     + h[idx_h_RRRR]  )  +
                                      (f[idx_f_L]      + f[idx_f_RRRRR] )*(g[idx_g_L]      + g[idx_g_RRRRR] )*(h[idx_h_L]      + h[idx_h_RRRRR] )) +
                    quarter*(d_coef_f*(f[idx_f_LLLLLL] + f[idx_f_R]     )*(g[idx_g_LLLLLL] + g[idx_g_R]     )*(h[idx_h_LLLLLL] + h[idx_h_R]     )  +
                                      (f[idx_f_LLLLL]  + f[idx_f_RR]    )*(g[idx_g_LLLLL]  + g[idx_g_RR]    )*(h[idx_h_LLLLL]  + h[idx_h_RR]    )  +
                                      (f[idx_f_LLLL]   + f[idx_f_RRR]   )*(g[idx_g_LLLL]   + g[idx_g_RRR]   )*(h[idx_h_LLLL]   + h[idx_h_RRR]   )  +
                                      (f[idx_f_LLL]    + f[idx_f_RRRR]  )*(g[idx_g_LLL]    + g[idx_g_RRRR]  )*(h[idx_h_LLL]    + h[idx_h_RRRR]  )  +
                                      (f[idx_f_LL]     + f[idx_f_RRRRR] )*(g[idx_g_LL]     + g[idx_g_RRRRR] )*(h[idx_h_LL]     + h[idx_h_RRRRR] )  +
                                      (f[idx_f_L]      + f[idx_f_RRRRRR])*(g[idx_g_L]      + g[idx_g_RRRRRR])*(h[idx_h_L]      + h[idx_h_RRRRRR])));
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        /*
         * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
         */
        
        const int num_ghosts_0_f = num_ghosts_f[0];
        const int num_ghosts_1_f = num_ghosts_f[1];
        const int ghostcell_dim_0_f = ghostcell_dims_f[0];
        
        const int num_ghosts_0_g = num_ghosts_g[0];
        const int num_ghosts_1_g = num_ghosts_g[1];
        const int ghostcell_dim_0_g = ghostcell_dims_g[0];
        
        const int num_ghosts_0_h = num_ghosts_h[0];
        const int num_ghosts_1_h = num_ghosts_h[1];
        const int ghostcell_dim_0_h = ghostcell_dims_h[0];
        
        if (d_stencil_width == 3)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_L = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_L = (i - 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_R = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_h_L = (i - 1 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_R = (i + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    F_face_x[idx_face_x] += dt*(
                        quarter*(d_coef_a*(f[idx_f_L] + f[idx_f_R])*(g[idx_g_L] + g[idx_g_R])*(h[idx_h_L] + h[idx_h_R])));
                }
            }
        }
        else if (d_stencil_width == 5)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LL = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L  = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R  = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_LL = (i - 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_L  = (i - 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_R  = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RR = (i + 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_h_LL = (i - 2 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_L  = (i - 1 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_R  = (i + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RR = (i + 1 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    F_face_x[idx_face_x] += dt*(
                        quarter*(d_coef_a*(f[idx_f_L]  + f[idx_f_R] )*(g[idx_g_L]  + g[idx_g_R] )*(h[idx_h_L]  + h[idx_h_R] )) +
                        quarter*(d_coef_b*(f[idx_f_LL] + f[idx_f_R] )*(g[idx_g_LL] + g[idx_g_R] )*(h[idx_h_LL] + h[idx_h_R] )  +
                                          (f[idx_f_L]  + f[idx_f_RR])*(g[idx_g_L]  + g[idx_g_RR])*(h[idx_h_L]  + h[idx_h_RR])));
                }
            }
        }
        else if (d_stencil_width == 7)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LLL = (i - 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LL  = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L   = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R   = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR  = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRR = (i + 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_LLL = (i - 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LL  = (i - 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_L   = (i - 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_R   = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RR  = (i + 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRR = (i + 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_h_LLL = (i - 3 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_LL  = (i - 2 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_L   = (i - 1 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_R   = (i + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RR  = (i + 1 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RRR = (i + 2 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    F_face_x[idx_face_x] += dt*(
                        quarter*(d_coef_a*(f[idx_f_L]   + f[idx_f_R]  )*(g[idx_g_L]   + g[idx_g_R]  )*(h[idx_h_L]   + h[idx_h_R]  )) +
                        quarter*(d_coef_b*(f[idx_f_LL]  + f[idx_f_R]  )*(g[idx_g_LL]  + g[idx_g_R]  )*(h[idx_h_LL]  + h[idx_h_R]  )  +
                                          (f[idx_f_L]   + f[idx_f_RR] )*(g[idx_g_L]   + g[idx_g_RR] )*(h[idx_h_L]   + h[idx_h_RR] )) +
                        quarter*(d_coef_c*(f[idx_f_LLL] + f[idx_f_R]  )*(g[idx_g_LLL] + g[idx_g_R]  )*(h[idx_h_LLL] + h[idx_h_R]  )  +
                                          (f[idx_f_LL]  + f[idx_f_RR] )*(g[idx_g_LL]  + g[idx_g_RR] )*(h[idx_h_LL]  + h[idx_h_RR] )  +
                                          (f[idx_f_L]   + f[idx_f_RRR])*(g[idx_g_L]   + g[idx_g_RRR])*(h[idx_h_L]   + h[idx_h_RRR])));
                }
            }
        }
        else if (d_stencil_width == 9)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LLLL = (i - 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLL  = (i - 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LL   = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L    = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R    = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR   = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRR  = (i + 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRR = (i + 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_LLLL = (i - 4 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LLL  = (i - 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LL   = (i - 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_L    = (i - 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_R    = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RR   = (i + 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRR  = (i + 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRRR = (i + 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_h_LLLL = (i - 4 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_LLL  = (i - 3 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_LL   = (i - 2 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_L    = (i - 1 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_R    = (i + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RR   = (i + 1 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RRR  = (i + 2 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RRRR = (i + 3 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    F_face_x[idx_face_x] += dt*(
                        quarter*(d_coef_a*(f[idx_f_L]    + f[idx_f_R]   )*(g[idx_g_L]    + g[idx_g_R]   )*(h[idx_h_L]    + h[idx_h_R]   )) +
                        quarter*(d_coef_b*(f[idx_f_LL]   + f[idx_f_R]   )*(g[idx_g_LL]   + g[idx_g_R]   )*(h[idx_h_LL]   + h[idx_h_R]   )  +
                                          (f[idx_f_L]    + f[idx_f_RR]  )*(g[idx_g_L]    + g[idx_g_RR]  )*(h[idx_h_L]    + h[idx_h_RR]  )) +
                        quarter*(d_coef_c*(f[idx_f_LLL]  + f[idx_f_R]   )*(g[idx_g_LLL]  + g[idx_g_R]   )*(h[idx_h_LLL]  + h[idx_h_R]   )  +
                                          (f[idx_f_LL]   + f[idx_f_RR]  )*(g[idx_g_LL]   + g[idx_g_RR]  )*(h[idx_h_LL]   + h[idx_h_RR]  )  +
                                          (f[idx_f_L]    + f[idx_f_RRR] )*(g[idx_g_L]    + g[idx_g_RRR] )*(h[idx_h_L]    + h[idx_h_RRR] )) +
                        quarter*(d_coef_d*(f[idx_f_LLLL] + f[idx_f_R]   )*(g[idx_g_LLLL] + g[idx_g_R]   )*(h[idx_h_LLLL] + h[idx_h_R]   )  +
                                          (f[idx_f_LLL]  + f[idx_f_RR]  )*(g[idx_g_LLL]  + g[idx_g_RR]  )*(h[idx_h_LLL]  + h[idx_h_RR]  )  +
                                          (f[idx_f_LL]   + f[idx_f_RRR] )*(g[idx_g_LL]   + g[idx_g_RRR] )*(h[idx_h_LL]   + h[idx_h_RRR] )  +
                                          (f[idx_f_L]    + f[idx_f_RRRR])*(g[idx_g_L]    + g[idx_g_RRRR])*(h[idx_h_L]    + h[idx_h_RRRR])));
                }
            }
        }
        else if (d_stencil_width == 11)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LLLLL = (i - 5 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLLL  = (i - 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLL   = (i - 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LL    = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L     = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R     = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR    = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRR   = (i + 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRR  = (i + 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRRR = (i + 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_LLLLL = (i - 5 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LLLL  = (i - 4 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LLL   = (i - 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LL    = (i - 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_L     = (i - 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_R     = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RR    = (i + 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRR   = (i + 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRRR  = (i + 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRRRR = (i + 4 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_h_LLLLL = (i - 5 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_LLLL  = (i - 4 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_LLL   = (i - 3 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_LL    = (i - 2 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_L     = (i - 1 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_R     = (i + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RR    = (i + 1 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RRR   = (i + 2 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RRRR  = (i + 3 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RRRRR = (i + 4 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    F_face_x[idx_face_x] += dt*(
                        quarter*(d_coef_a*(f[idx_f_L]     + f[idx_f_R]    )*(g[idx_g_L]     + g[idx_g_R]    )*(h[idx_h_L]     + h[idx_h_R]    )) +
                        quarter*(d_coef_b*(f[idx_f_LL]    + f[idx_f_R]    )*(g[idx_g_LL]    + g[idx_g_R]    )*(h[idx_h_LL]    + h[idx_h_R]    )  +
                                          (f[idx_f_L]     + f[idx_f_RR]   )*(g[idx_g_L]     + g[idx_g_RR]   )*(h[idx_h_L]     + h[idx_h_RR]   )) +
                        quarter*(d_coef_c*(f[idx_f_LLL]   + f[idx_f_R]    )*(g[idx_g_LLL]   + g[idx_g_R]    )*(h[idx_h_LLL]   + h[idx_h_R]    )  +
                                          (f[idx_f_LL]    + f[idx_f_RR]   )*(g[idx_g_LL]    + g[idx_g_RR]   )*(h[idx_h_LL]    + h[idx_h_RR]   )  +
                                          (f[idx_f_L]     + f[idx_f_RRR]  )*(g[idx_g_L]     + g[idx_g_RRR]  )*(h[idx_h_L]     + h[idx_h_RRR]  )) +
                        quarter*(d_coef_d*(f[idx_f_LLLL]  + f[idx_f_R]    )*(g[idx_g_LLLL]  + g[idx_g_R]    )*(h[idx_h_LLLL]  + h[idx_h_R]    )  +
                                          (f[idx_f_LLL]   + f[idx_f_RR]   )*(g[idx_g_LLL]   + g[idx_g_RR]   )*(h[idx_h_LLL]   + h[idx_h_RR]   )  +
                                          (f[idx_f_LL]    + f[idx_f_RRR]  )*(g[idx_g_LL]    + g[idx_g_RRR]  )*(h[idx_h_LL]    + h[idx_h_RRR]  )  +
                                          (f[idx_f_L]     + f[idx_f_RRRR] )*(g[idx_g_L]     + g[idx_g_RRRR] )*(h[idx_h_L]     + h[idx_h_RRRR] )) +
                        quarter*(d_coef_e*(f[idx_f_LLLLL] + f[idx_f_R]    )*(g[idx_g_LLLLL] + g[idx_g_R]    )*(h[idx_h_LLLLL] + h[idx_h_R]    )  +
                                          (f[idx_f_LLLL]  + f[idx_f_RR]   )*(g[idx_g_LLLL]  + g[idx_g_RR]   )*(h[idx_h_LLLL]  + h[idx_h_RR]   )  +
                                          (f[idx_f_LLL]   + f[idx_f_RRR]  )*(g[idx_g_LLL]   + g[idx_g_RRR]  )*(h[idx_h_LLL]   + h[idx_h_RRR]  )  +
                                          (f[idx_f_LL]    + f[idx_f_RRRR] )*(g[idx_g_LL]    + g[idx_g_RRRR] )*(h[idx_h_LL]    + h[idx_h_RRRR] )  +
                                          (f[idx_f_L]     + f[idx_f_RRRRR])*(g[idx_g_L]     + g[idx_g_RRRRR])*(h[idx_h_L]     + h[idx_h_RRRRR])));
                }
            }
        }
        else if (d_stencil_width == 13)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    const int idx_face_x = i +
                        j*(interior_dim_0 + 1);
                    
                    const int idx_f_LLLLLL = (i - 6 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLLLL  = (i - 5 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLLL   = (i - 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LLL    = (i - 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_LL     = (i - 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_L      = (i - 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_R      = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RR     = (i + 1 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRR    = (i + 2 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRR   = (i + 3 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRRR  = (i + 4 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_RRRRRR = (i + 5 + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_LLLLLL = (i - 6 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LLLLL  = (i - 5 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LLLL   = (i - 4 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LLL    = (i - 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_LL     = (i - 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_L      = (i - 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_R      = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RR     = (i + 1 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRR    = (i + 2 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRRR   = (i + 3 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRRRR  = (i + 4 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_RRRRRR = (i + 5 + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_h_LLLLLL = (i - 6 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_LLLLL  = (i - 5 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_LLLL   = (i - 4 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_LLL    = (i - 3 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_LL     = (i - 2 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_L      = (i - 1 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_R      = (i + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RR     = (i + 1 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RRR    = (i + 2 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RRRR   = (i + 3 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RRRRR  = (i + 4 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_RRRRRR = (i + 5 + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    F_face_x[idx_face_x] += dt*(
                        quarter*(d_coef_a*(f[idx_f_L]      + f[idx_f_R]     )*(g[idx_g_L]      + g[idx_g_R]     )*(h[idx_h_L]      + h[idx_h_R]     )) +
                        quarter*(d_coef_b*(f[idx_f_LL]     + f[idx_f_R]     )*(g[idx_g_LL]     + g[idx_g_R]     )*(h[idx_h_LL]     + h[idx_h_R]     )  +
                                          (f[idx_f_L]      + f[idx_f_RR]    )*(g[idx_g_L]      + g[idx_g_RR]    )*(h[idx_h_L]      + h[idx_h_RR]    )) +
                        quarter*(d_coef_c*(f[idx_f_LLL]    + f[idx_f_R]     )*(g[idx_g_LLL]    + g[idx_g_R]     )*(h[idx_h_LLL]    + h[idx_h_R]     )  +
                                          (f[idx_f_LL]     + f[idx_f_RR]    )*(g[idx_g_LL]     + g[idx_g_RR]    )*(h[idx_h_LL]     + h[idx_h_RR]    )  +
                                          (f[idx_f_L]      + f[idx_f_RRR]   )*(g[idx_g_L]      + g[idx_g_RRR]   )*(h[idx_h_L]      + h[idx_h_RRR]   )) +
                        quarter*(d_coef_d*(f[idx_f_LLLL]   + f[idx_f_R]     )*(g[idx_g_LLLL]   + g[idx_g_R]     )*(h[idx_h_LLLL]   + h[idx_h_R]     )  +
                                          (f[idx_f_LLL]    + f[idx_f_RR]    )*(g[idx_g_LLL]    + g[idx_g_RR]    )*(h[idx_h_LLL]    + h[idx_h_RR]    )  +
                                          (f[idx_f_LL]     + f[idx_f_RRR]   )*(g[idx_g_LL]     + g[idx_g_RRR]   )*(h[idx_h_LL]     + h[idx_h_RRR]   )  +
                                          (f[idx_f_L]      + f[idx_f_RRRR]  )*(g[idx_g_L]      + g[idx_g_RRRR]  )*(h[idx_h_L]      + h[idx_h_RRRR]  )) +
                        quarter*(d_coef_e*(f[idx_f_LLLLL]  + f[idx_f_R]     )*(g[idx_g_LLLLL]  + g[idx_g_R]     )*(h[idx_h_LLLLL]  + h[idx_h_R]     )  +
                                          (f[idx_f_LLLL]   + f[idx_f_RR]    )*(g[idx_g_LLLL]   + g[idx_g_RR]    )*(h[idx_h_LLLL]   + h[idx_h_RR]    )  +
                                          (f[idx_f_LLL]    + f[idx_f_RRR]   )*(g[idx_g_LLL]    + g[idx_g_RRR]   )*(h[idx_h_LLL]    + h[idx_h_RRR]   )  +
                                          (f[idx_f_LL]     + f[idx_f_RRRR]  )*(g[idx_g_LL]     + g[idx_g_RRRR]  )*(h[idx_h_LL]     + h[idx_h_RRRR]  )  +
                                          (f[idx_f_L]      + f[idx_f_RRRRR] )*(g[idx_g_L]      + g[idx_g_RRRRR] )*(h[idx_h_L]      + h[idx_h_RRRRR] )) +
                        quarter*(d_coef_f*(f[idx_f_LLLLLL] + f[idx_f_R]     )*(g[idx_g_LLLLLL] + g[idx_g_R]     )*(h[idx_h_LLLLLL] + h[idx_h_R]     )  +
                                          (f[idx_f_LLLLL]  + f[idx_f_RR]    )*(g[idx_g_LLLLL]  + g[idx_g_RR]    )*(h[idx_h_LLLLL]  + h[idx_h_RR]    )  +
                                          (f[idx_f_LLLL]   + f[idx_f_RRR]   )*(g[idx_g_LLLL]   + g[idx_g_RRR]   )*(h[idx_h_LLLL]   + h[idx_h_RRR]   )  +
                                          (f[idx_f_LLL]    + f[idx_f_RRRR]  )*(g[idx_g_LLL]    + g[idx_g_RRRR]  )*(h[idx_h_LLL]    + h[idx_h_RRRR]  )  +
                                          (f[idx_f_LL]     + f[idx_f_RRRRR] )*(g[idx_g_LL]     + g[idx_g_RRRRR] )*(h[idx_h_LL]     + h[idx_h_RRRRR] )  +
                                          (f[idx_f_L]      + f[idx_f_RRRRRR])*(g[idx_g_L]      + g[idx_g_RRRRRR])*(h[idx_h_L]      + h[idx_h_RRRRRR])));
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        /*
         * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
         */
        
        const int num_ghosts_0_f = num_ghosts_f[0];
        const int num_ghosts_1_f = num_ghosts_f[1];
        const int num_ghosts_2_f = num_ghosts_f[2];
        const int ghostcell_dim_0_f = ghostcell_dims_f[0];
        const int ghostcell_dim_1_f = ghostcell_dims_f[1];
        
        const int num_ghosts_0_g = num_ghosts_g[0];
        const int num_ghosts_1_g = num_ghosts_g[1];
        const int num_ghosts_2_g = num_ghosts_g[2];
        const int ghostcell_dim_0_g = ghostcell_dims_g[0];
        const int ghostcell_dim_1_g = ghostcell_dims_g[1];
        
        const int num_ghosts_0_h = num_ghosts_h[0];
        const int num_ghosts_1_h = num_ghosts_h[1];
        const int num_ghosts_2_h = num_ghosts_h[2];
        const int ghostcell_dim_0_h = ghostcell_dims_h[0];
        const int ghostcell_dim_1_h = ghostcell_dims_h[1];
        
        if (d_stencil_width == 3)
        {
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
                        
                        const int idx_f_L = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_L = (i - 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_R = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_L = (i - 1 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h + 
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_R = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_x[idx_face_x] += dt*(
                            quarter*(d_coef_a*(f[idx_f_L] + f[idx_f_R])*(g[idx_g_L] + g[idx_g_R])*(h[idx_h_L] + h[idx_h_R])));
                    }
                }
            }
        }
        else if (d_stencil_width == 5)
        {
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
                        
                        const int idx_f_LL = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L  = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_LL = (i - 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_L  = (i - 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_R  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RR = (i + 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_LL = (i - 2 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_L  = (i - 1 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h + 
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_R  = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RR = (i + 1 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_x[idx_face_x] += dt*(
                            quarter*(d_coef_a*(f[idx_f_L]  + f[idx_f_R] )*(g[idx_g_L]  + g[idx_g_R] )*(h[idx_h_L]  + h[idx_h_R] )) +
                            quarter*(d_coef_b*(f[idx_f_LL] + f[idx_f_R] )*(g[idx_g_LL] + g[idx_g_R] )*(h[idx_h_LL] + h[idx_h_R] )  +
                                              (f[idx_f_L]  + f[idx_f_RR])*(g[idx_g_L]  + g[idx_g_RR])*(h[idx_h_L]  + h[idx_h_RR])));
                    }
                }
            }
        }
        else if (d_stencil_width == 7)
        {
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
                        
                        const int idx_f_LLL = (i - 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LL  = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L   = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR  = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRR = (i + 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_LLL = (i - 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LL  = (i - 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_L   = (i - 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_R   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RR  = (i + 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRR = (i + 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_LLL = (i - 3 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_LL  = (i - 2 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_L   = (i - 1 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h + 
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_R   = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RR  = (i + 1 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RRR = (i + 2 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_x[idx_face_x] += dt*(
                            quarter*(d_coef_a*(f[idx_f_L]   + f[idx_f_R]  )*(g[idx_g_L]   + g[idx_g_R]  )*(h[idx_h_L]   + h[idx_h_R]  )) +
                            quarter*(d_coef_b*(f[idx_f_LL]  + f[idx_f_R]  )*(g[idx_g_LL]  + g[idx_g_R]  )*(h[idx_h_LL]  + h[idx_h_R]  )  +
                                              (f[idx_f_L]   + f[idx_f_RR] )*(g[idx_g_L]   + g[idx_g_RR] )*(h[idx_h_L]   + h[idx_h_RR] )) +
                            quarter*(d_coef_c*(f[idx_f_LLL] + f[idx_f_R]  )*(g[idx_g_LLL] + g[idx_g_R]  )*(h[idx_h_LLL] + h[idx_h_R]  )  +
                                              (f[idx_f_LL]  + f[idx_f_RR] )*(g[idx_g_LL]  + g[idx_g_RR] )*(h[idx_h_LL]  + h[idx_h_RR] )  +
                                              (f[idx_f_L]   + f[idx_f_RRR])*(g[idx_g_L]   + g[idx_g_RRR])*(h[idx_h_L]   + h[idx_h_RRR])));
                    }
                }
            }
        }
        else if (d_stencil_width == 9)
        {
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
                        
                        const int idx_f_LLLL = (i - 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLL  = (i - 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LL   = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L    = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR   = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRR  = (i + 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRR = (i + 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_LLLL = (i - 4 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LLL  = (i - 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LL   = (i - 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_L    = (i - 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_R    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RR   = (i + 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRR  = (i + 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRRR = (i + 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_LLLL = (i - 4 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h + 
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_LLL  = (i - 3 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_LL   = (i - 2 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_L    = (i - 1 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h + 
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_R    = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RR   = (i + 1 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RRR  = (i + 2 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RRRR = (i + 3 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_x[idx_face_x] += dt*(
                            quarter*(d_coef_a*(f[idx_f_L]    + f[idx_f_R]   )*(g[idx_g_L]    + g[idx_g_R]   )*(h[idx_h_L]    + h[idx_h_R]   )) +
                            quarter*(d_coef_b*(f[idx_f_LL]   + f[idx_f_R]   )*(g[idx_g_LL]   + g[idx_g_R]   )*(h[idx_h_LL]   + h[idx_h_R]   )  +
                                              (f[idx_f_L]    + f[idx_f_RR]  )*(g[idx_g_L]    + g[idx_g_RR]  )*(h[idx_h_L]    + h[idx_h_RR]  )) +
                            quarter*(d_coef_c*(f[idx_f_LLL]  + f[idx_f_R]   )*(g[idx_g_LLL]  + g[idx_g_R]   )*(h[idx_h_LLL]  + h[idx_h_R]   )  +
                                              (f[idx_f_LL]   + f[idx_f_RR]  )*(g[idx_g_LL]   + g[idx_g_RR]  )*(h[idx_h_LL]   + h[idx_h_RR]  )  +
                                              (f[idx_f_L]    + f[idx_f_RRR] )*(g[idx_g_L]    + g[idx_g_RRR] )*(h[idx_h_L]    + h[idx_h_RRR] )) +
                            quarter*(d_coef_d*(f[idx_f_LLLL] + f[idx_f_R]   )*(g[idx_g_LLLL] + g[idx_g_R]   )*(h[idx_h_LLLL] + h[idx_h_R]   )  +
                                              (f[idx_f_LLL]  + f[idx_f_RR]  )*(g[idx_g_LLL]  + g[idx_g_RR]  )*(h[idx_h_LLL]  + h[idx_h_RR]  )  +
                                              (f[idx_f_LL]   + f[idx_f_RRR] )*(g[idx_g_LL]   + g[idx_g_RRR] )*(h[idx_h_LL]   + h[idx_h_RRR] )  +
                                              (f[idx_f_L]    + f[idx_f_RRRR])*(g[idx_g_L]    + g[idx_g_RRRR])*(h[idx_h_L]    + h[idx_h_RRRR])));
                    }
                }
            }
        }
        else if (d_stencil_width == 11)
        {
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
                        
                        const int idx_f_LLLLL = (i - 5 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLLL  = (i - 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLL   = (i - 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LL    = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L     = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR    = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRR   = (i + 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRR  = (i + 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRRR = (i + 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_LLLLL = (i - 5 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LLLL  = (i - 4 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LLL   = (i - 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LL    = (i - 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_L     = (i - 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_R     = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RR    = (i + 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRR   = (i + 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRRR  = (i + 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRRRR = (i + 4 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_LLLLL = (i - 5 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h + 
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_LLLL  = (i - 4 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h + 
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_LLL   = (i - 3 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_LL    = (i - 2 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_L     = (i - 1 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h + 
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_R     = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RR    = (i + 1 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RRR   = (i + 2 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RRRR  = (i + 3 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RRRRR = (i + 4 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_x[idx_face_x] += dt*(
                            quarter*(d_coef_a*(f[idx_f_L]     + f[idx_f_R]    )*(g[idx_g_L]     + g[idx_g_R]    )*(h[idx_h_L]     + h[idx_h_R]    )) +
                            quarter*(d_coef_b*(f[idx_f_LL]    + f[idx_f_R]    )*(g[idx_g_LL]    + g[idx_g_R]    )*(h[idx_h_LL]    + h[idx_h_R]    )  +
                                              (f[idx_f_L]     + f[idx_f_RR]   )*(g[idx_g_L]     + g[idx_g_RR]   )*(h[idx_h_L]     + h[idx_h_RR]   )) +
                            quarter*(d_coef_c*(f[idx_f_LLL]   + f[idx_f_R]    )*(g[idx_g_LLL]   + g[idx_g_R]    )*(h[idx_h_LLL]   + h[idx_h_R]    )  +
                                              (f[idx_f_LL]    + f[idx_f_RR]   )*(g[idx_g_LL]    + g[idx_g_RR]   )*(h[idx_h_LL]    + h[idx_h_RR]   )  +
                                              (f[idx_f_L]     + f[idx_f_RRR]  )*(g[idx_g_L]     + g[idx_g_RRR]  )*(h[idx_h_L]     + h[idx_h_RRR]  )) +
                            quarter*(d_coef_d*(f[idx_f_LLLL]  + f[idx_f_R]    )*(g[idx_g_LLLL]  + g[idx_g_R]    )*(h[idx_h_LLLL]  + h[idx_h_R]    )  +
                                              (f[idx_f_LLL]   + f[idx_f_RR]   )*(g[idx_g_LLL]   + g[idx_g_RR]   )*(h[idx_h_LLL]   + h[idx_h_RR]   )  +
                                              (f[idx_f_LL]    + f[idx_f_RRR]  )*(g[idx_g_LL]    + g[idx_g_RRR]  )*(h[idx_h_LL]    + h[idx_h_RRR]  )  +
                                              (f[idx_f_L]     + f[idx_f_RRRR] )*(g[idx_g_L]     + g[idx_g_RRRR] )*(h[idx_h_L]     + h[idx_h_RRRR] )) +
                            quarter*(d_coef_e*(f[idx_f_LLLLL] + f[idx_f_R]    )*(g[idx_g_LLLLL] + g[idx_g_R]    )*(h[idx_h_LLLLL] + h[idx_h_R]    )  +
                                              (f[idx_f_LLLL]  + f[idx_f_RR]   )*(g[idx_g_LLLL]  + g[idx_g_RR]   )*(h[idx_h_LLLL]  + h[idx_h_RR]   )  +
                                              (f[idx_f_LLL]   + f[idx_f_RRR]  )*(g[idx_g_LLL]   + g[idx_g_RRR]  )*(h[idx_h_LLL]   + h[idx_h_RRR]  )  +
                                              (f[idx_f_LL]    + f[idx_f_RRRR] )*(g[idx_g_LL]    + g[idx_g_RRRR] )*(h[idx_h_LL]    + h[idx_h_RRRR] )  +
                                              (f[idx_f_L]     + f[idx_f_RRRRR])*(g[idx_g_L]     + g[idx_g_RRRRR])*(h[idx_h_L]     + h[idx_h_RRRRR])));
                    }
                }
            }
        }
        else if (d_stencil_width == 13)
        {
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
                        
                        const int idx_f_LLLLLL = (i - 6 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLLLL  = (i - 5 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLLL   = (i - 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LLL    = (i - 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_LL     = (i - 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_L      = (i - 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f + 
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_R      = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RR     = (i + 1 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRR    = (i + 2 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRR   = (i + 3 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRRR  = (i + 4 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_RRRRRR = (i + 5 + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_LLLLLL = (i - 6 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LLLLL  = (i - 5 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LLLL   = (i - 4 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LLL    = (i - 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_LL     = (i - 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_L      = (i - 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g + 
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_R      = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RR     = (i + 1 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRR    = (i + 2 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRRR   = (i + 3 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRRRR  = (i + 4 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_RRRRRR = (i + 5 + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_LLLLLL = (i - 6 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h + 
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_LLLLL  = (i - 5 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h + 
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_LLLL   = (i - 4 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h + 
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_LLL    = (i - 3 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_LL     = (i - 2 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_L      = (i - 1 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h + 
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_R      = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RR     = (i + 1 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RRR    = (i + 2 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RRRR   = (i + 3 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RRRRR  = (i + 4 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_RRRRRR = (i + 5 + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_x[idx_face_x] += dt*(
                            quarter*(d_coef_a*(f[idx_f_L]      + f[idx_f_R]     )*(g[idx_g_L]      + g[idx_g_R]     )*(h[idx_h_L]      + h[idx_h_R]     )) +
                            quarter*(d_coef_b*(f[idx_f_LL]     + f[idx_f_R]     )*(g[idx_g_LL]     + g[idx_g_R]     )*(h[idx_h_LL]     + h[idx_h_R]     )  +
                                              (f[idx_f_L]      + f[idx_f_RR]    )*(g[idx_g_L]      + g[idx_g_RR]    )*(h[idx_h_L]      + h[idx_h_RR]    )) +
                            quarter*(d_coef_c*(f[idx_f_LLL]    + f[idx_f_R]     )*(g[idx_g_LLL]    + g[idx_g_R]     )*(h[idx_h_LLL]    + h[idx_h_R]     )  +
                                              (f[idx_f_LL]     + f[idx_f_RR]    )*(g[idx_g_LL]     + g[idx_g_RR]    )*(h[idx_h_LL]     + h[idx_h_RR]    )  +
                                              (f[idx_f_L]      + f[idx_f_RRR]   )*(g[idx_g_L]      + g[idx_g_RRR]   )*(h[idx_h_L]      + h[idx_h_RRR]   )) +
                            quarter*(d_coef_d*(f[idx_f_LLLL]   + f[idx_f_R]     )*(g[idx_g_LLLL]   + g[idx_g_R]     )*(h[idx_h_LLLL]   + h[idx_h_R]     )  +
                                              (f[idx_f_LLL]    + f[idx_f_RR]    )*(g[idx_g_LLL]    + g[idx_g_RR]    )*(h[idx_h_LLL]    + h[idx_h_RR]    )  +
                                              (f[idx_f_LL]     + f[idx_f_RRR]   )*(g[idx_g_LL]     + g[idx_g_RRR]   )*(h[idx_h_LL]     + h[idx_h_RRR]   )  +
                                              (f[idx_f_L]      + f[idx_f_RRRR]  )*(g[idx_g_L]      + g[idx_g_RRRR]  )*(h[idx_h_L]      + h[idx_h_RRRR]  )) +
                            quarter*(d_coef_e*(f[idx_f_LLLLL]  + f[idx_f_R]     )*(g[idx_g_LLLLL]  + g[idx_g_R]     )*(h[idx_h_LLLLL]  + h[idx_h_R]     )  +
                                              (f[idx_f_LLLL]   + f[idx_f_RR]    )*(g[idx_g_LLLL]   + g[idx_g_RR]    )*(h[idx_h_LLLL]   + h[idx_h_RR]    )  +
                                              (f[idx_f_LLL]    + f[idx_f_RRR]   )*(g[idx_g_LLL]    + g[idx_g_RRR]   )*(h[idx_h_LLL]    + h[idx_h_RRR]   )  +
                                              (f[idx_f_LL]     + f[idx_f_RRRR]  )*(g[idx_g_LL]     + g[idx_g_RRRR]  )*(h[idx_h_LL]     + h[idx_h_RRRR]  )  +
                                              (f[idx_f_L]      + f[idx_f_RRRRR] )*(g[idx_g_L]      + g[idx_g_RRRRR] )*(h[idx_h_L]      + h[idx_h_RRRRR] )) +
                            quarter*(d_coef_f*(f[idx_f_LLLLLL] + f[idx_f_R]     )*(g[idx_g_LLLLLL] + g[idx_g_R]     )*(h[idx_h_LLLLLL] + h[idx_h_R]     )  +
                                              (f[idx_f_LLLLL]  + f[idx_f_RR]    )*(g[idx_g_LLLLL]  + g[idx_g_RR]    )*(h[idx_h_LLLLL]  + h[idx_h_RR]    )  +
                                              (f[idx_f_LLLL]   + f[idx_f_RRR]   )*(g[idx_g_LLLL]   + g[idx_g_RRR]   )*(h[idx_h_LLLL]   + h[idx_h_RRR]   )  +
                                              (f[idx_f_LLL]    + f[idx_f_RRRR]  )*(g[idx_g_LLL]    + g[idx_g_RRRR]  )*(h[idx_h_LLL]    + h[idx_h_RRRR]  )  +
                                              (f[idx_f_LL]     + f[idx_f_RRRRR] )*(g[idx_g_LL]     + g[idx_g_RRRRR] )*(h[idx_h_LL]     + h[idx_h_RRRRR] )  +
                                              (f[idx_f_L]      + f[idx_f_RRRRRR])*(g[idx_g_L]      + g[idx_g_RRRRRR])*(h[idx_h_L]      + h[idx_h_RRRRRR])));
                    }
                }
            }
        }
    }
}
