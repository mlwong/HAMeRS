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
    d_use_DRP4 = d_convective_flux_reconstructor_db->
        getBoolWithDefault("use_DRP4", false);
    d_use_DRP4 = d_convective_flux_reconstructor_db->
        getBoolWithDefault("d_use_DRP4", d_use_DRP4);
    
    if (d_use_DRP4)
    {
        d_stencil_width = d_convective_flux_reconstructor_db->
            getIntegerWithDefault("stencil_width", 9);
        d_stencil_width = d_convective_flux_reconstructor_db->
            getIntegerWithDefault("d_stencil_width", d_stencil_width);
        d_order = 4;
        
        if (d_stencil_width < 9 || d_stencil_width > 13)
        {
            TBOX_ERROR("ConvectiveFluxReconstructorKEP::computeConvectiveFluxAndSourceOnPatch:"
                " Only 9-point, 11-point, 13-point stencil KEP schemes are implemented!");
        }
    }
    else
    {
        d_order = d_convective_flux_reconstructor_db->
            getIntegerWithDefault("order", 2);
        d_order = d_convective_flux_reconstructor_db->
            getIntegerWithDefault("d_order", d_order);
        
        if (d_order == 2)
        {
            d_stencil_width = 3;
        }
        else if (d_order == 4)
        {
            d_stencil_width = 5;
        }
        else if (d_order == 6)
        {
            d_stencil_width = 7;
        }
        else if (d_order == 8)
        {
            d_stencil_width = 9;
        }
        else if (d_order == 10)
        {
            d_stencil_width = 11;
        }
        else if (d_order == 12)
        {
            d_stencil_width = 13;
        }
        else
        {
            TBOX_ERROR("ConvectiveFluxReconstructorKEP::computeConvectiveFluxAndSourceOnPatch:"
                " Only 3-point, 5-point, 7-point, 9-point, 11-point, 13-point stencil central schemes are implemented!");
        }
    }
    
    d_coef_a = double(0);
    d_coef_b = double(0);
    d_coef_c = double(0);
    d_coef_d = double(0);
    d_coef_e = double(0);
    d_coef_f = double(0);
    
    if (d_use_DRP4)
    {
        if (d_stencil_width == 9)
        {
            d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*4;
            
            d_coef_a = double( 0.841570125482);
            d_coef_b = double(-0.244678631765);
            d_coef_c = double( 0.059463584768);
            d_coef_d = double(-0.007650904064);
        }
        else if (d_stencil_width == 11)
        {
            d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*5;
            
            d_coef_a = double( 0.872756993962);
            d_coef_b = double(-0.286511173973);
            d_coef_c = double( 0.090320001280);
            d_coef_d = double(-0.020779405824);
            d_coef_e = double( 0.002484594688);
        }
        else if (d_stencil_width == 13)
        {
            d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*6;
            
            d_coef_a = double( 0.907646591371);
            d_coef_b = double(-0.337048393268);
            d_coef_c = double( 0.133442885327);
            d_coef_d = double(-0.045246480208);
            d_coef_e = double( 0.011169294114);
            d_coef_f = double(-0.001456501759);
        }
    }
    else
    {
        // Use traditional central schemes.
        
        if (d_stencil_width == 3)
        {
            d_num_conv_ghosts = hier::IntVector::getOne(d_dim);
            
            d_coef_a = double(1)/double(2);
        }
        else if (d_stencil_width == 5)
        {
            d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*2;
            
            d_coef_a =  double(2)/double(3);
            d_coef_b = -double(1)/double(12);
        }
        else if (d_stencil_width == 7)
        {
            d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*3;
            
            d_coef_a =  double(3)/double(4);
            d_coef_b = -double(3)/double(20);
            d_coef_c =  double(1)/double(60);
        }
        else if (d_stencil_width == 9)
        {
            d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*4;
            
            d_coef_a =  double(4)/double(5);
            d_coef_b = -double(1)/double(5);
            d_coef_c =  double(4)/double(105);
            d_coef_d = -double(1)/double(280);
        }
        else if (d_stencil_width == 11)
        {
            d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*5;
            
            d_coef_a =  double(5)/double(6);
            d_coef_b = -double(5)/double(21);
            d_coef_c =  double(5)/double(84);
            d_coef_d = -double(5)/double(504);
            d_coef_e =  double(1)/double(1260);
        }
        else if (d_stencil_width == 13)
        {
            d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*6;
            
            d_coef_a =  double(6)/double(7);
            d_coef_b = -double(15)/double(56);
            d_coef_c =  double(5)/double(63);
            d_coef_d = -double(1)/double(56);
            d_coef_e =  double(3)/double(1155);
            d_coef_f = -double(1)/double(5544);
        }
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
    restart_db->putBool("d_use_DRP4", d_use_DRP4);
    restart_db->putInteger("d_stencil_width", d_stencil_width);
    restart_db->putInteger("d_order", d_order);
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
    if (d_flow_model_type != FLOW_MODEL::SINGLE_SPECIES &&
        d_flow_model_type != FLOW_MODEL::FOUR_EQN_CONSERVATIVE &&
        d_flow_model_type != FLOW_MODEL::FIVE_EQN_ALLAIRE)
    {
        TBOX_ERROR("ConvectiveFluxReconstructorKEP::computeConvectiveFluxAndSourceOnPatch:"
            " KEP schemes can only be used for flow models: SINGLE_SPECIES, FOUR_EQN_CONSERVATIVE or FIVE_EQN_ALLAIRE!");
    }
    
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
    
    // Initialize the flux to be zero.
    
    convective_flux->fill(double(0));
    
    /*
     * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
     */
    
    d_flow_model->registerPatchWithDataContext(patch, data_context);
    
    std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
    
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("DENSITY", d_num_conv_ghosts));
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRESSURE", d_num_conv_ghosts));
    num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("TOTAL_ENERGY", d_num_conv_ghosts));
    
    if (d_flow_model_type == FLOW_MODEL::FOUR_EQN_CONSERVATIVE || d_flow_model_type == FLOW_MODEL::FIVE_EQN_ALLAIRE)
    {
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", d_num_conv_ghosts));
    }
    
    if (d_flow_model_type == FLOW_MODEL::FIVE_EQN_ALLAIRE)
    {
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VOLUME_FRACTIONS", d_num_conv_ghosts));
    }
    
    d_flow_model->registerDerivedVariables(num_subghosts_of_data);
    
    d_flow_model->allocateMemoryForDerivedCellData();
    
    d_flow_model->computeDerivedCellData();
    
    /*
     * Get the pointers to the velocity and convective flux cell data inside the flow model.
     * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
     */
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > density      = d_flow_model->getCellData("DENSITY");
    HAMERS_SHARED_PTR<pdat::CellData<double> > velocity     = d_flow_model->getCellData("VELOCITY");
    HAMERS_SHARED_PTR<pdat::CellData<double> > pressure     = d_flow_model->getCellData("PRESSURE");
    HAMERS_SHARED_PTR<pdat::CellData<double> > total_energy = d_flow_model->getCellData("TOTAL_ENERGY");
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > mass_fractions;
    HAMERS_SHARED_PTR<pdat::CellData<double> > volume_fractions;
    
    if (d_flow_model_type == FLOW_MODEL::FOUR_EQN_CONSERVATIVE || d_flow_model_type == FLOW_MODEL::FIVE_EQN_ALLAIRE)
    {
        mass_fractions = d_flow_model->getCellData("MASS_FRACTIONS");
    }
    
    if (d_flow_model_type == FLOW_MODEL::FIVE_EQN_ALLAIRE)
    {
        volume_fractions = d_flow_model->getCellData("VOLUME_FRACTIONS");
    }
    
    hier::IntVector num_subghosts_density = density->getGhostCellWidth();
    hier::IntVector subghostcell_dims_density = density->getGhostBox().numberCells();
    
    hier::IntVector num_subghosts_pressure = pressure->getGhostCellWidth();
    hier::IntVector subghostcell_dims_pressure = pressure->getGhostBox().numberCells();
    
    hier::IntVector num_subghosts_total_energy = total_energy->getGhostCellWidth();
    hier::IntVector subghostcell_dims_total_energy = total_energy->getGhostBox().numberCells();
    
    // Initialize data container for specific total enthalpy.
    
    HAMERS_SHARED_PTR<pdat::CellData<double> > specific_total_enthalpy(new pdat::CellData<double>(
       interior_box, 1, d_num_conv_ghosts));
    
    // Get the pointers to the cell data.
    
    double* rho = density->getPointer(0);
    double* p   = pressure->getPointer(0);
    double* E   = total_energy->getPointer(0);
    double* H   = specific_total_enthalpy->getPointer(0);
    
    // Compute the specific total enthalpy.
    
    t_reconstruct_flux->start();
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the interior dimension.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_conv_ghosts_0 = d_num_conv_ghosts[0];
        
        const int num_subghosts_0_density = num_subghosts_density[0];
        const int num_subghosts_0_pressure = num_subghosts_pressure[0];
        const int num_subghosts_0_total_energy = num_subghosts_total_energy[0];
        
#ifdef HAMERS_ENABLE_SIMD
        #pragma omp simd
#endif
        for (int i = -num_conv_ghosts_0; i < interior_dim_0 + num_conv_ghosts_0; i++)
        {
            // Compute the linear indices.
            
            const int idx = i + num_conv_ghosts_0;
            
            const int idx_density      = i + num_subghosts_0_density;
            const int idx_pressure     = i + num_subghosts_0_pressure;
            const int idx_total_energy = i + num_subghosts_0_total_energy;
            
            H[idx] = (E[idx_total_energy] + p[idx_pressure])/(rho[idx_density]);
        }
        
        
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the interior dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        /*
         * Get the numbers of ghost cells and the dimensions of the ghost cell boxes.
         */
        
        const int num_conv_ghosts_0 = d_num_conv_ghosts[0];
        const int num_conv_ghosts_1 = d_num_conv_ghosts[1];
        
        const int num_subghosts_0_density = num_subghosts_density[0];
        const int num_subghosts_1_density = num_subghosts_density[1];
        const int subghostcell_dim_0_density = subghostcell_dims_density[0];
        
        const int num_subghosts_0_pressure = num_subghosts_pressure[0];
        const int num_subghosts_1_pressure = num_subghosts_pressure[1];
        const int subghostcell_dim_0_pressure = subghostcell_dims_pressure[0];
        
        const int num_subghosts_0_total_energy = num_subghosts_total_energy[0];
        const int num_subghosts_1_total_energy = num_subghosts_total_energy[1];
        const int subghostcell_dim_0_total_energy = subghostcell_dims_total_energy[0];
        
        for (int j = -num_conv_ghosts_1; j < interior_dim_1 + num_conv_ghosts_1; j++)
        {
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd
#endif
            for (int i = -num_conv_ghosts_0; i < interior_dim_0 + num_conv_ghosts_0; i++)
            {
                // Compute the linear indices.
                
                const int idx = (i + num_conv_ghosts_0) +
                    (j + num_conv_ghosts_1)*(interior_dim_0 + 2*num_conv_ghosts_0);
                
                const int idx_density = (i + num_subghosts_0_density) +
                    (j + num_subghosts_1_density)*subghostcell_dim_0_density;
                
                const int idx_pressure = (i + num_subghosts_0_pressure) +
                    (j + num_subghosts_1_pressure)*subghostcell_dim_0_pressure;
                
                const int idx_total_energy = (i + num_subghosts_0_total_energy) +
                    (j + num_subghosts_1_total_energy)*subghostcell_dim_0_total_energy;
                
                H[idx] = (E[idx_total_energy] + p[idx_pressure])/(rho[idx_density]);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the interior dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_conv_ghosts_0 = d_num_conv_ghosts[0];
        const int num_conv_ghosts_1 = d_num_conv_ghosts[1];
        const int num_conv_ghosts_2 = d_num_conv_ghosts[2];
        
        const int num_subghosts_0_density = num_subghosts_density[0];
        const int num_subghosts_1_density = num_subghosts_density[1];
        const int num_subghosts_2_density = num_subghosts_density[2];
        const int subghostcell_dim_0_density = subghostcell_dims_density[0];
        const int subghostcell_dim_1_density = subghostcell_dims_density[1];
        
        const int num_subghosts_0_pressure = num_subghosts_pressure[0];
        const int num_subghosts_1_pressure = num_subghosts_pressure[1];
        const int num_subghosts_2_pressure = num_subghosts_pressure[2];
        const int subghostcell_dim_0_pressure = subghostcell_dims_pressure[0];
        const int subghostcell_dim_1_pressure = subghostcell_dims_pressure[1];
        
        const int num_subghosts_0_total_energy = num_subghosts_total_energy[0];
        const int num_subghosts_1_total_energy = num_subghosts_total_energy[1];
        const int num_subghosts_2_total_energy = num_subghosts_total_energy[2];
        const int subghostcell_dim_0_total_energy = subghostcell_dims_total_energy[0];
        const int subghostcell_dim_1_total_energy = subghostcell_dims_total_energy[1];
        
        for (int k = -num_conv_ghosts_2; k < interior_dim_2 + num_conv_ghosts_2; k++)
        {
            for (int j = -num_conv_ghosts_1; j < interior_dim_1 + num_conv_ghosts_1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = -num_conv_ghosts_0; i < interior_dim_0 + num_conv_ghosts_0; i++)
                {
                    // Compute the linear indices.
                    
                    const int idx = (i + num_conv_ghosts_0) +
                        (j + num_conv_ghosts_1)*(interior_dim_0 + 2*num_conv_ghosts_0) +
                        (k + num_conv_ghosts_2)*(interior_dim_0 + 2*num_conv_ghosts_0)*
                            (interior_dim_1 + 2*num_conv_ghosts_1);
                    
                    const int idx_density = (i + num_subghosts_0_density) +
                        (j + num_subghosts_1_density)*subghostcell_dim_0_density +
                        (k + num_subghosts_2_density)*subghostcell_dim_0_density*
                            subghostcell_dim_1_density;
                    
                    const int idx_pressure = (i + num_subghosts_0_pressure) +
                        (j + num_subghosts_1_pressure)*subghostcell_dim_0_pressure +
                        (k + num_subghosts_2_pressure)*subghostcell_dim_0_pressure*
                            subghostcell_dim_1_pressure;
                    
                    const int idx_total_energy = (i + num_subghosts_0_total_energy) +
                        (j + num_subghosts_1_total_energy)*subghostcell_dim_0_total_energy +
                        (k + num_subghosts_2_total_energy)*subghostcell_dim_0_total_energy*
                            subghostcell_dim_1_total_energy;
                    
                    H[idx] = (E[idx_total_energy] + p[idx_pressure])/(rho[idx_density]);
                }
            }
        }
    }
    
    t_reconstruct_flux->stop();
    
    /*
     * Compute the flux.
     */
    
    t_reconstruct_flux->start();
    
    if (d_dim == tbox::Dimension(1))
    {
        if (d_flow_model_type == FLOW_MODEL::SINGLE_SPECIES)
        {
            /*
             * Reconstruct flux in x-direction.
             */
            
            // Density equation.
            
            addQuadraticTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                0,
                0,
                0,
                dt);
            
            // Momentum equation.
            
            addLinearTermToConvectiveFluxX(convective_flux,
                pressure,
                1,
                0,
                dt);
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                velocity,
                1,
                0,
                0,
                0,
                dt);
            
            // Total energy equation.
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                specific_total_enthalpy,
                2,
                0,
                0,
                0,
                dt);
        }
        else if (d_flow_model_type == FLOW_MODEL::FOUR_EQN_CONSERVATIVE ||
                 d_flow_model_type == FLOW_MODEL::FIVE_EQN_ALLAIRE)
        {
            const int num_species = d_flow_model->getNumberOfSpecies();
            
            /*
             * Reconstruct flux in x-direction.
             */
            
            // Partial density equations.
            
            for (int si = 0; si < num_species; si++)
            {
                addCubicTermToConvectiveFluxX(convective_flux,
                    density,
                    velocity,
                    mass_fractions,
                    si,
                    0,
                    0,
                    si,
                    dt);
            }
            
            // Momentum equation.
            
            addLinearTermToConvectiveFluxX(convective_flux,
                pressure,
                num_species,
                0,
                dt);
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                velocity,
                num_species,
                0,
                0,
                0,
                dt);
            
            // Total energy equation.
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                specific_total_enthalpy,
                num_species + 1,
                0,
                0,
                0,
                dt);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        if (d_flow_model_type == FLOW_MODEL::SINGLE_SPECIES)
        {
            /*
             * Reconstruct flux in x-direction.
             */
            
            // Density equation.
            
            addQuadraticTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                0,
                0,
                0,
                dt);
            
            // Momentum equation in x-direction.
            
            addLinearTermToConvectiveFluxX(convective_flux,
                pressure,
                1,
                0,
                dt);
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                velocity,
                1,
                0,
                0,
                0,
                dt);
            
            // Momentum equation in y-direction.
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                velocity,
                2,
                0,
                1,
                0,
                dt);
            
            // Total energy equation.
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                specific_total_enthalpy,
                3,
                0,
                0,
                0,
                dt);
            
            /*
             * Reconstruct flux in y-direction.
             */
            
            // Density equation.
            
            addQuadraticTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                0,
                0,
                1,
                dt);
            
            // Momentum equation in x-direction.
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                velocity,
                1,
                0,
                0,
                1,
                dt);
            
            // Momentum equation in y-direction.
            
            addLinearTermToConvectiveFluxY(convective_flux,
                pressure,
                2,
                0,
                dt);
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                velocity,
                2,
                0,
                1,
                1,
                dt);
            
            // Total energy equation.
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                specific_total_enthalpy,
                3,
                0,
                1,
                0,
                dt);
        }
        else if (d_flow_model_type == FLOW_MODEL::FOUR_EQN_CONSERVATIVE ||
                 d_flow_model_type == FLOW_MODEL::FIVE_EQN_ALLAIRE)
        {
            const int num_species = d_flow_model->getNumberOfSpecies();
            
            /*
             * Reconstruct flux in x-direction.
             */
            
            // Partial density equations.
            
            for (int si = 0; si < num_species; si++)
            {
                addCubicTermToConvectiveFluxX(convective_flux,
                    density,
                    velocity,
                    mass_fractions,
                    si,
                    0,
                    0,
                    si,
                    dt);
            }
            
            // Momentum equation in x-direction.
            
            addLinearTermToConvectiveFluxX(convective_flux,
                pressure,
                num_species,
                0,
                dt);
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                velocity,
                num_species,
                0,
                0,
                0,
                dt);
            
            // Momentum equation in y-direction.
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                velocity,
                num_species + 1,
                0,
                1,
                0,
                dt);
            
            // Total energy equation.
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                specific_total_enthalpy,
                num_species + 2,
                0,
                0,
                0,
                dt);
            
            /*
             * Reconstruct flux in y-direction.
             */
            
            // Partial density equations.
            
            for (int si = 0; si < num_species; si++)
            {
                addCubicTermToConvectiveFluxY(convective_flux,
                    density,
                    velocity,
                    mass_fractions,
                    si,
                    0,
                    1,
                    si,
                    dt);
            }
            
            // Momentum equation in x-direction.
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                velocity,
                num_species,
                0,
                0,
                1,
                dt);
            
            // Momentum equation in y-direction.
            
            addLinearTermToConvectiveFluxY(convective_flux,
                pressure,
                num_species + 1,
                0,
                dt);
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                velocity,
                num_species + 1,
                0,
                1,
                1,
                dt);
            
            // Total energy equation.
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                specific_total_enthalpy,
                num_species + 2,
                0,
                1,
                0,
                dt);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        if (d_flow_model_type == FLOW_MODEL::SINGLE_SPECIES)
        {
            /*
             * Reconstruct flux in x-direction.
             */
            
            // Density equation.
            
            addQuadraticTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                0,
                0,
                0,
                dt);
            
            // Momentum equation in x-direction.
            
            addLinearTermToConvectiveFluxX(convective_flux,
                pressure,
                1,
                0,
                dt);
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                velocity,
                1,
                0,
                0,
                0,
                dt);
            
            // Momentum equation in y-direction.
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                velocity,
                2,
                0,
                1,
                0,
                dt);
            
            // Momentum equation in z-direction.
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                velocity,
                3,
                0,
                2,
                0,
                dt);
            
            // Total energy equation.
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                specific_total_enthalpy,
                4,
                0,
                0,
                0,
                dt);
            
            /*
             * Reconstruct flux in y-direction.
             */
            
            // Density equation.
            
            addQuadraticTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                0,
                0,
                1,
                dt);
            
            // Momentum equation in x-direction.
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                velocity,
                1,
                0,
                0,
                1,
                dt);
            
            // Momentum equation in y-direction.
            
            addLinearTermToConvectiveFluxY(convective_flux,
                pressure,
                2,
                0,
                dt);
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                velocity,
                2,
                0,
                1,
                1,
                dt);
            
            // Momentum equation in z-direction.
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                velocity,
                3,
                0,
                2,
                1,
                dt);
            
            // Total energy equation.
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                specific_total_enthalpy,
                4,
                0,
                1,
                0,
                dt);
            
            /*
             * Reconstruct flux in z-direction.
             */
            
            // Density equation.
            
            addQuadraticTermToConvectiveFluxZ(convective_flux,
                density,
                velocity,
                0,
                0,
                2,
                dt);
            
            // Momentum equation in x-direction.
            
            addCubicTermToConvectiveFluxZ(convective_flux,
                density,
                velocity,
                velocity,
                1,
                0,
                0,
                2,
                dt);
            
            // Momentum equation in y-direction.
            
            addCubicTermToConvectiveFluxZ(convective_flux,
                density,
                velocity,
                velocity,
                2,
                0,
                1,
                2,
                dt);
            
            // Momentum equation in z-direction.
            
            addLinearTermToConvectiveFluxZ(convective_flux,
                pressure,
                3,
                0,
                dt);
            
            addCubicTermToConvectiveFluxZ(convective_flux,
                density,
                velocity,
                velocity,
                3,
                0,
                2,
                2,
                dt);
            
            // Total energy equation.
            
            addCubicTermToConvectiveFluxZ(convective_flux,
                density,
                velocity,
                specific_total_enthalpy,
                4,
                0,
                2,
                0,
                dt);
        }
        else if (d_flow_model_type == FLOW_MODEL::FOUR_EQN_CONSERVATIVE ||
                 d_flow_model_type == FLOW_MODEL::FIVE_EQN_ALLAIRE)
        {
            const int num_species = d_flow_model->getNumberOfSpecies();
            
            /*
             * Reconstruct flux in x-direction.
             */
             
            // Partial density equations.
            
            for (int si = 0; si < num_species; si++)
            {
                addCubicTermToConvectiveFluxX(convective_flux,
                    density,
                    velocity,
                    mass_fractions,
                    si,
                    0,
                    0,
                    si,
                    dt);
            }
            
            // Momentum equation in x-direction.
            
            addLinearTermToConvectiveFluxX(convective_flux,
                pressure,
                num_species,
                0,
                dt);
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                velocity,
                num_species,
                0,
                0,
                0,
                dt);
            
            // Momentum equation in y-direction.
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                velocity,
                num_species + 1,
                0,
                1,
                0,
                dt);
            
            // Momentum equation in z-direction.
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                velocity,
                num_species + 2,
                0,
                2,
                0,
                dt);
            
            // Total energy equation.
            
            addCubicTermToConvectiveFluxX(convective_flux,
                density,
                velocity,
                specific_total_enthalpy,
                num_species + 3,
                0,
                0,
                0,
                dt);
            
            /*
             * Reconstruct flux in y-direction.
             */
            
            // Partial density equations.
            
            for (int si = 0; si < num_species; si++)
            {
                addCubicTermToConvectiveFluxY(convective_flux,
                    density,
                    velocity,
                    mass_fractions,
                    si,
                    0,
                    1,
                    si,
                    dt);
            }
            
            // Momentum equation in x-direction.
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                velocity,
                num_species,
                0,
                0,
                1,
                dt);
            
            // Momentum equation in y-direction.
            
            addLinearTermToConvectiveFluxY(convective_flux,
                pressure,
                num_species + 1,
                0,
                dt);
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                velocity,
                num_species + 1,
                0,
                1,
                1,
                dt);
            
            // Momentum equation in z-direction.
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                velocity,
                num_species + 2,
                0,
                2,
                1,
                dt);
            
            // Total energy equation.
            
            addCubicTermToConvectiveFluxY(convective_flux,
                density,
                velocity,
                specific_total_enthalpy,
                num_species + 3,
                0,
                1,
                0,
                dt);
            
            /*
             * Reconstruct flux in z-direction.
             */
            
            // Partial density equations.
            
            for (int si = 0; si < num_species; si++)
            {
                addCubicTermToConvectiveFluxZ(convective_flux,
                    density,
                    velocity,
                    mass_fractions,
                    si,
                    0,
                    2,
                    si,
                    dt);
            }
            
            // Momentum equation in x-direction.
            
            addCubicTermToConvectiveFluxZ(convective_flux,
                density,
                velocity,
                velocity,
                num_species,
                0,
                0,
                2,
                dt);
            
            // Momentum equation in y-direction.
            
            addCubicTermToConvectiveFluxZ(convective_flux,
                density,
                velocity,
                velocity,
                num_species + 1,
                0,
                1,
                2,
                dt);
            
            // Momentum equation in z-direction.
            
            addLinearTermToConvectiveFluxZ(convective_flux,
                pressure,
                num_species + 2,
                0,
                dt);
            
            addCubicTermToConvectiveFluxZ(convective_flux,
                density,
                velocity,
                velocity,
                num_species + 2,
                0,
                2,
                2,
                dt);
            
            // Total energy equation.
            
            addCubicTermToConvectiveFluxZ(convective_flux,
                density,
                velocity,
                specific_total_enthalpy,
                num_species + 3,
                0,
                2,
                0,
                dt);
        }
    }
    
    t_reconstruct_flux->stop();
    
    /*
     * Compute the source.
     */
    
    t_compute_source->start();
    
    if (d_flow_model_type == FLOW_MODEL::FIVE_EQN_ALLAIRE)
    {
        addSourceTermsToVolumeFractionEquations(
            source,
            velocity,
            volume_fractions,
            dx,
            dt);
    }
    
    t_compute_source->stop();
    
    /*
     * Unregister the patch and data of all registered derived cell variables in the flow model.
     */
    
    d_flow_model->unregisterPatch();
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
         * Get the dimension.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Get the numbers of ghost cells.
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
                    d_coef_a*(f[idx_f_L] + f[idx_f_R]));
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
                    d_coef_a*((f[idx_f_L]  + f[idx_f_R] )) +
                    d_coef_b*((f[idx_f_LL] + f[idx_f_R] )  +
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
                    d_coef_a*((f[idx_f_L]   + f[idx_f_R]  )) +
                    d_coef_b*((f[idx_f_LL]  + f[idx_f_R]  )  +
                              (f[idx_f_L]   + f[idx_f_RR] )) +
                    d_coef_c*((f[idx_f_LLL] + f[idx_f_R]  )  +
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
                    d_coef_a*((f[idx_f_L]    + f[idx_f_R]   )) +
                    d_coef_b*((f[idx_f_LL]   + f[idx_f_R]   )  +
                              (f[idx_f_L]    + f[idx_f_RR]  )) +
                    d_coef_c*((f[idx_f_LLL]  + f[idx_f_R]   )  +
                              (f[idx_f_LL]   + f[idx_f_RR]  )  +
                              (f[idx_f_L]    + f[idx_f_RRR] )) +
                    d_coef_d*((f[idx_f_LLLL] + f[idx_f_R]   )  +
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
                    d_coef_a*((f[idx_f_L]     + f[idx_f_R]    )) +
                    d_coef_b*((f[idx_f_LL]    + f[idx_f_R]    )  +
                              (f[idx_f_L]     + f[idx_f_RR]   )) +
                    d_coef_c*((f[idx_f_LLL]   + f[idx_f_R]    )  +
                              (f[idx_f_LL]    + f[idx_f_RR]   )  +
                              (f[idx_f_L]     + f[idx_f_RRR]  )) +
                    d_coef_d*((f[idx_f_LLLL]  + f[idx_f_R]    )  +
                              (f[idx_f_LLL]   + f[idx_f_RR]   )  +
                              (f[idx_f_LL]    + f[idx_f_RRR]  )  +
                              (f[idx_f_L]     + f[idx_f_RRRR] )) +
                    d_coef_e*((f[idx_f_LLLLL] + f[idx_f_R]    )  +
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
                    d_coef_a*((f[idx_f_L]      + f[idx_f_R]     )) +
                    d_coef_b*((f[idx_f_LL]     + f[idx_f_R]     )  +
                              (f[idx_f_L]      + f[idx_f_RR]    )) +
                    d_coef_c*((f[idx_f_LLL]    + f[idx_f_R]     )  +
                              (f[idx_f_LL]     + f[idx_f_RR]    )  +
                              (f[idx_f_L]      + f[idx_f_RRR]   )) +
                    d_coef_d*((f[idx_f_LLLL]   + f[idx_f_R]     )  +
                              (f[idx_f_LLL]    + f[idx_f_RR]    )  +
                              (f[idx_f_LL]     + f[idx_f_RRR]   )  +
                              (f[idx_f_L]      + f[idx_f_RRRR]  )) +
                    d_coef_e*((f[idx_f_LLLLL]  + f[idx_f_R]     )  +
                              (f[idx_f_LLLL]   + f[idx_f_RR]    )  +
                              (f[idx_f_LLL]    + f[idx_f_RRR]   )  +
                              (f[idx_f_LL]     + f[idx_f_RRRR]  )  +
                              (f[idx_f_L]      + f[idx_f_RRRRR] )) +
                    d_coef_f*((f[idx_f_LLLLLL] + f[idx_f_R]     )  +
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
         * Get the numbers of ghost cells and the dimensions of the ghost cell box.
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
                        d_coef_a*(f[idx_f_L] + f[idx_f_R]));
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
                        d_coef_a*((f[idx_f_L]  + f[idx_f_R] )) +
                        d_coef_b*((f[idx_f_LL] + f[idx_f_R] )  +
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
                        d_coef_a*((f[idx_f_L]   + f[idx_f_R]  )) +
                        d_coef_b*((f[idx_f_LL]  + f[idx_f_R]  )  +
                                  (f[idx_f_L]   + f[idx_f_RR] )) +
                        d_coef_c*((f[idx_f_LLL] + f[idx_f_R]  )  +
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
                        d_coef_a*((f[idx_f_L]    + f[idx_f_R]   )) +
                        d_coef_b*((f[idx_f_LL]   + f[idx_f_R]   )  +
                                  (f[idx_f_L]    + f[idx_f_RR]  )) +
                        d_coef_c*((f[idx_f_LLL]  + f[idx_f_R]   )  +
                                  (f[idx_f_LL]   + f[idx_f_RR]  )  +
                                  (f[idx_f_L]    + f[idx_f_RRR] )) +
                        d_coef_d*((f[idx_f_LLLL] + f[idx_f_R]   )  +
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
                        d_coef_a*((f[idx_f_L]     + f[idx_f_R]    )) +
                        d_coef_b*((f[idx_f_LL]    + f[idx_f_R]    )  +
                                  (f[idx_f_L]     + f[idx_f_RR]   )) +
                        d_coef_c*((f[idx_f_LLL]   + f[idx_f_R]    )  +
                                  (f[idx_f_LL]    + f[idx_f_RR]   )  +
                                  (f[idx_f_L]     + f[idx_f_RRR]  )) +
                        d_coef_d*((f[idx_f_LLLL]  + f[idx_f_R]    )  +
                                  (f[idx_f_LLL]   + f[idx_f_RR]   )  +
                                  (f[idx_f_LL]    + f[idx_f_RRR]  )  +
                                  (f[idx_f_L]     + f[idx_f_RRRR] )) +
                        d_coef_e*((f[idx_f_LLLLL] + f[idx_f_R]    )  +
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
                        d_coef_a*((f[idx_f_L]      + f[idx_f_R]     )) +
                        d_coef_b*((f[idx_f_LL]     + f[idx_f_R]     )  +
                                  (f[idx_f_L]      + f[idx_f_RR]    )) +
                        d_coef_c*((f[idx_f_LLL]    + f[idx_f_R]     )  +
                                  (f[idx_f_LL]     + f[idx_f_RR]    )  +
                                  (f[idx_f_L]      + f[idx_f_RRR]   )) +
                        d_coef_d*((f[idx_f_LLLL]   + f[idx_f_R]     )  +
                                  (f[idx_f_LLL]    + f[idx_f_RR]    )  +
                                  (f[idx_f_LL]     + f[idx_f_RRR]   )  +
                                  (f[idx_f_L]      + f[idx_f_RRRR]  )) +
                        d_coef_e*((f[idx_f_LLLLL]  + f[idx_f_R]     )  +
                                  (f[idx_f_LLLL]   + f[idx_f_RR]    )  +
                                  (f[idx_f_LLL]    + f[idx_f_RRR]   )  +
                                  (f[idx_f_LL]     + f[idx_f_RRRR]  )  +
                                  (f[idx_f_L]      + f[idx_f_RRRRR] )) +
                        d_coef_f*((f[idx_f_LLLLLL] + f[idx_f_R]     )  +
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
         * Get the numbers of ghost cells and the dimensions of the ghost cell box.
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
                            d_coef_a*(f[idx_f_L] + f[idx_f_R]));
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
                            d_coef_a*((f[idx_f_L]  + f[idx_f_R] )) +
                            d_coef_b*((f[idx_f_LL] + f[idx_f_R] )  +
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
                            d_coef_a*((f[idx_f_L]   + f[idx_f_R]  )) +
                            d_coef_b*((f[idx_f_LL]  + f[idx_f_R]  )  +
                                      (f[idx_f_L]   + f[idx_f_RR] )) +
                            d_coef_c*((f[idx_f_LLL] + f[idx_f_R]  )  +
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
                            d_coef_a*((f[idx_f_L]    + f[idx_f_R]   )) +
                            d_coef_b*((f[idx_f_LL]   + f[idx_f_R]   )  +
                                      (f[idx_f_L]    + f[idx_f_RR]  )) +
                            d_coef_c*((f[idx_f_LLL]  + f[idx_f_R]   )  +
                                      (f[idx_f_LL]   + f[idx_f_RR]  )  +
                                      (f[idx_f_L]    + f[idx_f_RRR] )) +
                            d_coef_d*((f[idx_f_LLLL] + f[idx_f_R]   )  +
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
                            d_coef_a*((f[idx_f_L]     + f[idx_f_R]    )) +
                            d_coef_b*((f[idx_f_LL]    + f[idx_f_R]    )  +
                                      (f[idx_f_L]     + f[idx_f_RR]   )) +
                            d_coef_c*((f[idx_f_LLL]   + f[idx_f_R]    )  +
                                      (f[idx_f_LL]    + f[idx_f_RR]   )  +
                                      (f[idx_f_L]     + f[idx_f_RRR]  )) +
                            d_coef_d*((f[idx_f_LLLL]  + f[idx_f_R]    )  +
                                      (f[idx_f_LLL]   + f[idx_f_RR]   )  +
                                      (f[idx_f_LL]    + f[idx_f_RRR]  )  +
                                      (f[idx_f_L]     + f[idx_f_RRRR] )) +
                            d_coef_e*((f[idx_f_LLLLL] + f[idx_f_R]    )  +
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
                            d_coef_a*((f[idx_f_L]      + f[idx_f_R]     )) +
                            d_coef_b*((f[idx_f_LL]     + f[idx_f_R]     )  +
                                      (f[idx_f_L]      + f[idx_f_RR]    )) +
                            d_coef_c*((f[idx_f_LLL]    + f[idx_f_R]     )  +
                                      (f[idx_f_LL]     + f[idx_f_RR]    )  +
                                      (f[idx_f_L]      + f[idx_f_RRR]   )) +
                            d_coef_d*((f[idx_f_LLLL]   + f[idx_f_R]     )  +
                                      (f[idx_f_LLL]    + f[idx_f_RR]    )  +
                                      (f[idx_f_LL]     + f[idx_f_RRR]   )  +
                                      (f[idx_f_L]      + f[idx_f_RRRR]  )) +
                            d_coef_e*((f[idx_f_LLLLL]  + f[idx_f_R]     )  +
                                      (f[idx_f_LLLL]   + f[idx_f_RR]    )  +
                                      (f[idx_f_LLL]    + f[idx_f_RRR]   )  +
                                      (f[idx_f_LL]     + f[idx_f_RRRR]  )  +
                                      (f[idx_f_L]      + f[idx_f_RRRRR] )) +
                            d_coef_f*((f[idx_f_LLLLLL] + f[idx_f_R]     )  +
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
         * Get the dimension.
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
                
                F_face_x[idx_face_x] += dt*half*(
                    d_coef_a*(f[idx_f_L] + f[idx_f_R])*(g[idx_g_L] + g[idx_g_R]));
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
                
                F_face_x[idx_face_x] += dt*half*(
                    d_coef_a*((f[idx_f_L]  + f[idx_f_R] )*(g[idx_g_L]  + g[idx_g_R] )) +
                    d_coef_b*((f[idx_f_LL] + f[idx_f_R] )*(g[idx_g_LL] + g[idx_g_R] )  +
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
                
                F_face_x[idx_face_x] += dt*half*(
                    d_coef_a*((f[idx_f_L]   + f[idx_f_R]  )*(g[idx_g_L]   + g[idx_g_R]  )) +
                    d_coef_b*((f[idx_f_LL]  + f[idx_f_R]  )*(g[idx_g_LL]  + g[idx_g_R]  )  +
                              (f[idx_f_L]   + f[idx_f_RR] )*(g[idx_g_L]   + g[idx_g_RR] )) +
                    d_coef_c*((f[idx_f_LLL] + f[idx_f_R]  )*(g[idx_g_LLL] + g[idx_g_R]  )  +
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
                
                F_face_x[idx_face_x] += dt*half*(
                    d_coef_a*((f[idx_f_L]    + f[idx_f_R]   )*(g[idx_g_L]    + g[idx_g_R]   )) +
                    d_coef_b*((f[idx_f_LL]   + f[idx_f_R]   )*(g[idx_g_LL]   + g[idx_g_R]   )  +
                              (f[idx_f_L]    + f[idx_f_RR]  )*(g[idx_g_L]    + g[idx_g_RR]  )) +
                    d_coef_c*((f[idx_f_LLL]  + f[idx_f_R]   )*(g[idx_g_LLL]  + g[idx_g_R]   )  +
                              (f[idx_f_LL]   + f[idx_f_RR]  )*(g[idx_g_LL]   + g[idx_g_RR]  )  +
                              (f[idx_f_L]    + f[idx_f_RRR] )*(g[idx_g_L]    + g[idx_g_RRR] )) +
                    d_coef_d*((f[idx_f_LLLL] + f[idx_f_R]   )*(g[idx_g_LLLL] + g[idx_g_R]   )  +
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
                
                F_face_x[idx_face_x] += dt*half*(
                    d_coef_a*((f[idx_f_L]     + f[idx_f_R]    )*(g[idx_g_L]     + g[idx_g_R]    )) +
                    d_coef_b*((f[idx_f_LL]    + f[idx_f_R]    )*(g[idx_g_LL]    + g[idx_g_R]    )  +
                              (f[idx_f_L]     + f[idx_f_RR]   )*(g[idx_g_L]     + g[idx_g_RR]   )) +
                    d_coef_c*((f[idx_f_LLL]   + f[idx_f_R]    )*(g[idx_g_LLL]   + g[idx_g_R]    )  +
                              (f[idx_f_LL]    + f[idx_f_RR]   )*(g[idx_g_LL]    + g[idx_g_RR]   )  +
                              (f[idx_f_L]     + f[idx_f_RRR]  )*(g[idx_g_L]     + g[idx_g_RRR]  )) +
                    d_coef_d*((f[idx_f_LLLL]  + f[idx_f_R]    )*(g[idx_g_LLLL]  + g[idx_g_R]    )  +
                              (f[idx_f_LLL]   + f[idx_f_RR]   )*(g[idx_g_LLL]   + g[idx_g_RR]   )  +
                              (f[idx_f_LL]    + f[idx_f_RRR]  )*(g[idx_g_LL]    + g[idx_g_RRR]  )  +
                              (f[idx_f_L]     + f[idx_f_RRRR] )*(g[idx_g_L]     + g[idx_g_RRRR] )) +
                    d_coef_e*((f[idx_f_LLLLL] + f[idx_f_R]    )*(g[idx_g_LLLLL] + g[idx_g_R]    )  +
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
                
                F_face_x[idx_face_x] += dt*half*(
                    d_coef_a*((f[idx_f_L]      + f[idx_f_R]     )*(g[idx_g_L]      + g[idx_g_R]     )) +
                    d_coef_b*((f[idx_f_LL]     + f[idx_f_R]     )*(g[idx_g_LL]     + g[idx_g_R]     )  +
                              (f[idx_f_L]      + f[idx_f_RR]    )*(g[idx_g_L]      + g[idx_g_RR]    )) +
                    d_coef_c*((f[idx_f_LLL]    + f[idx_f_R]     )*(g[idx_g_LLL]    + g[idx_g_R]     )  +
                              (f[idx_f_LL]     + f[idx_f_RR]    )*(g[idx_g_LL]     + g[idx_g_RR]    )  +
                              (f[idx_f_L]      + f[idx_f_RRR]   )*(g[idx_g_L]      + g[idx_g_RRR]   )) +
                    d_coef_d*((f[idx_f_LLLL]   + f[idx_f_R]     )*(g[idx_g_LLLL]   + g[idx_g_R]     )  +
                              (f[idx_f_LLL]    + f[idx_f_RR]    )*(g[idx_g_LLL]    + g[idx_g_RR]    )  +
                              (f[idx_f_LL]     + f[idx_f_RRR]   )*(g[idx_g_LL]     + g[idx_g_RRR]   )  +
                              (f[idx_f_L]      + f[idx_f_RRRR]  )*(g[idx_g_L]      + g[idx_g_RRRR]  )) +
                    d_coef_e*((f[idx_f_LLLLL]  + f[idx_f_R]     )*(g[idx_g_LLLLL]  + g[idx_g_R]     )  +
                              (f[idx_f_LLLL]   + f[idx_f_RR]    )*(g[idx_g_LLLL]   + g[idx_g_RR]    )  +
                              (f[idx_f_LLL]    + f[idx_f_RRR]   )*(g[idx_g_LLL]    + g[idx_g_RRR]   )  +
                              (f[idx_f_LL]     + f[idx_f_RRRR]  )*(g[idx_g_LL]     + g[idx_g_RRRR]  )  +
                              (f[idx_f_L]      + f[idx_f_RRRRR] )*(g[idx_g_L]      + g[idx_g_RRRRR] )) +
                    d_coef_f*((f[idx_f_LLLLLL] + f[idx_f_R]     )*(g[idx_g_LLLLLL] + g[idx_g_R]     )  +
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
                    
                    F_face_x[idx_face_x] += dt*half*(
                        d_coef_a*(f[idx_f_L] + f[idx_f_R])*(g[idx_g_L] + g[idx_g_R]));
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
                    
                    F_face_x[idx_face_x] += dt*half*(
                        d_coef_a*((f[idx_f_L]  + f[idx_f_R] )*(g[idx_g_L]  + g[idx_g_R] )) +
                        d_coef_b*((f[idx_f_LL] + f[idx_f_R] )*(g[idx_g_LL] + g[idx_g_R] )  +
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
                    
                    F_face_x[idx_face_x] += dt*half*(
                        d_coef_a*((f[idx_f_L]   + f[idx_f_R]  )*(g[idx_g_L]   + g[idx_g_R]  )) +
                        d_coef_b*((f[idx_f_LL]  + f[idx_f_R]  )*(g[idx_g_LL]  + g[idx_g_R]  )  +
                                  (f[idx_f_L]   + f[idx_f_RR] )*(g[idx_g_L]   + g[idx_g_RR] )) +
                        d_coef_c*((f[idx_f_LLL] + f[idx_f_R]  )*(g[idx_g_LLL] + g[idx_g_R]  )  +
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
                    
                    F_face_x[idx_face_x] += dt*half*(
                        d_coef_a*((f[idx_f_L]    + f[idx_f_R]   )*(g[idx_g_L]    + g[idx_g_R]   )) +
                        d_coef_b*((f[idx_f_LL]   + f[idx_f_R]   )*(g[idx_g_LL]   + g[idx_g_R]   )  +
                                  (f[idx_f_L]    + f[idx_f_RR]  )*(g[idx_g_L]    + g[idx_g_RR]  )) +
                        d_coef_c*((f[idx_f_LLL]  + f[idx_f_R]   )*(g[idx_g_LLL]  + g[idx_g_R]   )  +
                                  (f[idx_f_LL]   + f[idx_f_RR]  )*(g[idx_g_LL]   + g[idx_g_RR]  )  +
                                  (f[idx_f_L]    + f[idx_f_RRR] )*(g[idx_g_L]    + g[idx_g_RRR] )) +
                        d_coef_d*((f[idx_f_LLLL] + f[idx_f_R]   )*(g[idx_g_LLLL] + g[idx_g_R]   )  +
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
                    
                    F_face_x[idx_face_x] += dt*half*(
                        d_coef_a*((f[idx_f_L]     + f[idx_f_R]    )*(g[idx_g_L]     + g[idx_g_R]    )) +
                        d_coef_b*((f[idx_f_LL]    + f[idx_f_R]    )*(g[idx_g_LL]    + g[idx_g_R]    )  +
                                  (f[idx_f_L]     + f[idx_f_RR]   )*(g[idx_g_L]     + g[idx_g_RR]   )) +
                        d_coef_c*((f[idx_f_LLL]   + f[idx_f_R]    )*(g[idx_g_LLL]   + g[idx_g_R]    )  +
                                  (f[idx_f_LL]    + f[idx_f_RR]   )*(g[idx_g_LL]    + g[idx_g_RR]   )  +
                                  (f[idx_f_L]     + f[idx_f_RRR]  )*(g[idx_g_L]     + g[idx_g_RRR]  )) +
                        d_coef_d*((f[idx_f_LLLL]  + f[idx_f_R]    )*(g[idx_g_LLLL]  + g[idx_g_R]    )  +
                                  (f[idx_f_LLL]   + f[idx_f_RR]   )*(g[idx_g_LLL]   + g[idx_g_RR]   )  +
                                  (f[idx_f_LL]    + f[idx_f_RRR]  )*(g[idx_g_LL]    + g[idx_g_RRR]  )  +
                                  (f[idx_f_L]     + f[idx_f_RRRR] )*(g[idx_g_L]     + g[idx_g_RRRR] )) +
                        d_coef_e*((f[idx_f_LLLLL] + f[idx_f_R]    )*(g[idx_g_LLLLL] + g[idx_g_R]    )  +
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
                    
                    F_face_x[idx_face_x] += dt*half*(
                        d_coef_a*((f[idx_f_L]      + f[idx_f_R]     )*(g[idx_g_L]      + g[idx_g_R]     )) +
                        d_coef_b*((f[idx_f_LL]     + f[idx_f_R]     )*(g[idx_g_LL]     + g[idx_g_R]     )  +
                                  (f[idx_f_L]      + f[idx_f_RR]    )*(g[idx_g_L]      + g[idx_g_RR]    )) +
                        d_coef_c*((f[idx_f_LLL]    + f[idx_f_R]     )*(g[idx_g_LLL]    + g[idx_g_R]     )  +
                                  (f[idx_f_LL]     + f[idx_f_RR]    )*(g[idx_g_LL]     + g[idx_g_RR]    )  +
                                  (f[idx_f_L]      + f[idx_f_RRR]   )*(g[idx_g_L]      + g[idx_g_RRR]   )) +
                        d_coef_d*((f[idx_f_LLLL]   + f[idx_f_R]     )*(g[idx_g_LLLL]   + g[idx_g_R]     )  +
                                  (f[idx_f_LLL]    + f[idx_f_RR]    )*(g[idx_g_LLL]    + g[idx_g_RR]    )  +
                                  (f[idx_f_LL]     + f[idx_f_RRR]   )*(g[idx_g_LL]     + g[idx_g_RRR]   )  +
                                  (f[idx_f_L]      + f[idx_f_RRRR]  )*(g[idx_g_L]      + g[idx_g_RRRR]  )) +
                        d_coef_e*((f[idx_f_LLLLL]  + f[idx_f_R]     )*(g[idx_g_LLLLL]  + g[idx_g_R]     )  +
                                  (f[idx_f_LLLL]   + f[idx_f_RR]    )*(g[idx_g_LLLL]   + g[idx_g_RR]    )  +
                                  (f[idx_f_LLL]    + f[idx_f_RRR]   )*(g[idx_g_LLL]    + g[idx_g_RRR]   )  +
                                  (f[idx_f_LL]     + f[idx_f_RRRR]  )*(g[idx_g_LL]     + g[idx_g_RRRR]  )  +
                                  (f[idx_f_L]      + f[idx_f_RRRRR] )*(g[idx_g_L]      + g[idx_g_RRRRR] )) +
                        d_coef_f*((f[idx_f_LLLLLL] + f[idx_f_R]     )*(g[idx_g_LLLLLL] + g[idx_g_R]     )  +
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
                        
                        F_face_x[idx_face_x] += dt*half*(
                            d_coef_a*(f[idx_f_L] + f[idx_f_R])*(g[idx_g_L] + g[idx_g_R]));
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
                        
                        F_face_x[idx_face_x] += dt*half*(
                            d_coef_a*((f[idx_f_L]  + f[idx_f_R] )*(g[idx_g_L]  + g[idx_g_R] )) +
                            d_coef_b*((f[idx_f_LL] + f[idx_f_R] )*(g[idx_g_LL] + g[idx_g_R] )  +
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
                        
                        F_face_x[idx_face_x] += dt*half*(
                            d_coef_a*((f[idx_f_L]   + f[idx_f_R]  )*(g[idx_g_L]   + g[idx_g_R]  )) +
                            d_coef_b*((f[idx_f_LL]  + f[idx_f_R]  )*(g[idx_g_LL]  + g[idx_g_R]  )  +
                                      (f[idx_f_L]   + f[idx_f_RR] )*(g[idx_g_L]   + g[idx_g_RR] )) +
                            d_coef_c*((f[idx_f_LLL] + f[idx_f_R]  )*(g[idx_g_LLL] + g[idx_g_R]  )  +
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
                        
                        F_face_x[idx_face_x] += dt*half*(
                            d_coef_a*((f[idx_f_L]    + f[idx_f_R]   )*(g[idx_g_L]    + g[idx_g_R]   )) +
                            d_coef_b*((f[idx_f_LL]   + f[idx_f_R]   )*(g[idx_g_LL]   + g[idx_g_R]   )  +
                                      (f[idx_f_L]    + f[idx_f_RR]  )*(g[idx_g_L]    + g[idx_g_RR]  )) +
                            d_coef_c*((f[idx_f_LLL]  + f[idx_f_R]   )*(g[idx_g_LLL]  + g[idx_g_R]   )  +
                                      (f[idx_f_LL]   + f[idx_f_RR]  )*(g[idx_g_LL]   + g[idx_g_RR]  )  +
                                      (f[idx_f_L]    + f[idx_f_RRR] )*(g[idx_g_L]    + g[idx_g_RRR] )) +
                            d_coef_d*((f[idx_f_LLLL] + f[idx_f_R]   )*(g[idx_g_LLLL] + g[idx_g_R]   )  +
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
                        
                        F_face_x[idx_face_x] += dt*half*(
                            d_coef_a*((f[idx_f_L]     + f[idx_f_R]    )*(g[idx_g_L]     + g[idx_g_R]    )) +
                            d_coef_b*((f[idx_f_LL]    + f[idx_f_R]    )*(g[idx_g_LL]    + g[idx_g_R]    )  +
                                      (f[idx_f_L]     + f[idx_f_RR]   )*(g[idx_g_L]     + g[idx_g_RR]   )) +
                            d_coef_c*((f[idx_f_LLL]   + f[idx_f_R]    )*(g[idx_g_LLL]   + g[idx_g_R]    )  +
                                      (f[idx_f_LL]    + f[idx_f_RR]   )*(g[idx_g_LL]    + g[idx_g_RR]   )  +
                                      (f[idx_f_L]     + f[idx_f_RRR]  )*(g[idx_g_L]     + g[idx_g_RRR]  )) +
                            d_coef_d*((f[idx_f_LLLL]  + f[idx_f_R]    )*(g[idx_g_LLLL]  + g[idx_g_R]    )  +
                                      (f[idx_f_LLL]   + f[idx_f_RR]   )*(g[idx_g_LLL]   + g[idx_g_RR]   )  +
                                      (f[idx_f_LL]    + f[idx_f_RRR]  )*(g[idx_g_LL]    + g[idx_g_RRR]  )  +
                                      (f[idx_f_L]     + f[idx_f_RRRR] )*(g[idx_g_L]     + g[idx_g_RRRR] )) +
                            d_coef_e*((f[idx_f_LLLLL] + f[idx_f_R]    )*(g[idx_g_LLLLL] + g[idx_g_R]    )  +
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
                        
                        F_face_x[idx_face_x] += dt*half*(
                            d_coef_a*((f[idx_f_L]      + f[idx_f_R]     )*(g[idx_g_L]      + g[idx_g_R]     )) +
                            d_coef_b*((f[idx_f_LL]     + f[idx_f_R]     )*(g[idx_g_LL]     + g[idx_g_R]     )  +
                                      (f[idx_f_L]      + f[idx_f_RR]    )*(g[idx_g_L]      + g[idx_g_RR]    )) +
                            d_coef_c*((f[idx_f_LLL]    + f[idx_f_R]     )*(g[idx_g_LLL]    + g[idx_g_R]     )  +
                                      (f[idx_f_LL]     + f[idx_f_RR]    )*(g[idx_g_LL]     + g[idx_g_RR]    )  +
                                      (f[idx_f_L]      + f[idx_f_RRR]   )*(g[idx_g_L]      + g[idx_g_RRR]   )) +
                            d_coef_d*((f[idx_f_LLLL]   + f[idx_f_R]     )*(g[idx_g_LLLL]   + g[idx_g_R]     )  +
                                      (f[idx_f_LLL]    + f[idx_f_RR]    )*(g[idx_g_LLL]    + g[idx_g_RR]    )  +
                                      (f[idx_f_LL]     + f[idx_f_RRR]   )*(g[idx_g_LL]     + g[idx_g_RRR]   )  +
                                      (f[idx_f_L]      + f[idx_f_RRRR]  )*(g[idx_g_L]      + g[idx_g_RRRR]  )) +
                            d_coef_e*((f[idx_f_LLLLL]  + f[idx_f_R]     )*(g[idx_g_LLLLL]  + g[idx_g_R]     )  +
                                      (f[idx_f_LLLL]   + f[idx_f_RR]    )*(g[idx_g_LLLL]   + g[idx_g_RR]    )  +
                                      (f[idx_f_LLL]    + f[idx_f_RRR]   )*(g[idx_g_LLL]    + g[idx_g_RRR]   )  +
                                      (f[idx_f_LL]     + f[idx_f_RRRR]  )*(g[idx_g_LL]     + g[idx_g_RRRR]  )  +
                                      (f[idx_f_L]      + f[idx_f_RRRRR] )*(g[idx_g_L]      + g[idx_g_RRRRR] )) +
                            d_coef_f*((f[idx_f_LLLLLL] + f[idx_f_R]     )*(g[idx_g_LLLLLL] + g[idx_g_R]     )  +
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
         * Get the dimension.
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
                
                F_face_x[idx_face_x] += dt*quarter*(
                    d_coef_a*(f[idx_f_L] + f[idx_f_R])*(g[idx_g_L] + g[idx_g_R])*(h[idx_h_L] + h[idx_h_R]));
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
                
                F_face_x[idx_face_x] += dt*quarter*(
                    d_coef_a*((f[idx_f_L]  + f[idx_f_R] )*(g[idx_g_L]  + g[idx_g_R] )*(h[idx_h_L]  + h[idx_h_R] )) +
                    d_coef_b*((f[idx_f_LL] + f[idx_f_R] )*(g[idx_g_LL] + g[idx_g_R] )*(h[idx_h_LL] + h[idx_h_R] )  +
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
                
                F_face_x[idx_face_x] += dt*quarter*(
                    d_coef_a*((f[idx_f_L]   + f[idx_f_R]  )*(g[idx_g_L]   + g[idx_g_R]  )*(h[idx_h_L]   + h[idx_h_R]  )) +
                    d_coef_b*((f[idx_f_LL]  + f[idx_f_R]  )*(g[idx_g_LL]  + g[idx_g_R]  )*(h[idx_h_LL]  + h[idx_h_R]  )  +
                              (f[idx_f_L]   + f[idx_f_RR] )*(g[idx_g_L]   + g[idx_g_RR] )*(h[idx_h_L]   + h[idx_h_RR] )) +
                    d_coef_c*((f[idx_f_LLL] + f[idx_f_R]  )*(g[idx_g_LLL] + g[idx_g_R]  )*(h[idx_h_LLL] + h[idx_h_R]  )  +
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
                
                F_face_x[idx_face_x] += dt*quarter*(
                    d_coef_a*((f[idx_f_L]    + f[idx_f_R]   )*(g[idx_g_L]    + g[idx_g_R]   )*(h[idx_h_L]    + h[idx_h_R]   )) +
                    d_coef_b*((f[idx_f_LL]   + f[idx_f_R]   )*(g[idx_g_LL]   + g[idx_g_R]   )*(h[idx_h_LL]   + h[idx_h_R]   )  +
                              (f[idx_f_L]    + f[idx_f_RR]  )*(g[idx_g_L]    + g[idx_g_RR]  )*(h[idx_h_L]    + h[idx_h_RR]  )) +
                    d_coef_c*((f[idx_f_LLL]  + f[idx_f_R]   )*(g[idx_g_LLL]  + g[idx_g_R]   )*(h[idx_h_LLL]  + h[idx_h_R]   )  +
                              (f[idx_f_LL]   + f[idx_f_RR]  )*(g[idx_g_LL]   + g[idx_g_RR]  )*(h[idx_h_LL]   + h[idx_h_RR]  )  +
                              (f[idx_f_L]    + f[idx_f_RRR] )*(g[idx_g_L]    + g[idx_g_RRR] )*(h[idx_h_L]    + h[idx_h_RRR] )) +
                    d_coef_d*((f[idx_f_LLLL] + f[idx_f_R]   )*(g[idx_g_LLLL] + g[idx_g_R]   )*(h[idx_h_LLLL] + h[idx_h_R]   )  +
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
                
                F_face_x[idx_face_x] += dt*quarter*(
                    d_coef_a*((f[idx_f_L]     + f[idx_f_R]    )*(g[idx_g_L]     + g[idx_g_R]    )*(h[idx_h_L]     + h[idx_h_R]    )) +
                    d_coef_b*((f[idx_f_LL]    + f[idx_f_R]    )*(g[idx_g_LL]    + g[idx_g_R]    )*(h[idx_h_LL]    + h[idx_h_R]    )  +
                              (f[idx_f_L]     + f[idx_f_RR]   )*(g[idx_g_L]     + g[idx_g_RR]   )*(h[idx_h_L]     + h[idx_h_RR]   )) +
                    d_coef_c*((f[idx_f_LLL]   + f[idx_f_R]    )*(g[idx_g_LLL]   + g[idx_g_R]    )*(h[idx_h_LLL]   + h[idx_h_R]    )  +
                              (f[idx_f_LL]    + f[idx_f_RR]   )*(g[idx_g_LL]    + g[idx_g_RR]   )*(h[idx_h_LL]    + h[idx_h_RR]   )  +
                              (f[idx_f_L]     + f[idx_f_RRR]  )*(g[idx_g_L]     + g[idx_g_RRR]  )*(h[idx_h_L]     + h[idx_h_RRR]  )) +
                    d_coef_d*((f[idx_f_LLLL]  + f[idx_f_R]    )*(g[idx_g_LLLL]  + g[idx_g_R]    )*(h[idx_h_LLLL]  + h[idx_h_R]    )  +
                              (f[idx_f_LLL]   + f[idx_f_RR]   )*(g[idx_g_LLL]   + g[idx_g_RR]   )*(h[idx_h_LLL]   + h[idx_h_RR]   )  +
                              (f[idx_f_LL]    + f[idx_f_RRR]  )*(g[idx_g_LL]    + g[idx_g_RRR]  )*(h[idx_h_LL]    + h[idx_h_RRR]  )  +
                              (f[idx_f_L]     + f[idx_f_RRRR] )*(g[idx_g_L]     + g[idx_g_RRRR] )*(h[idx_h_L]     + h[idx_h_RRRR] )) +
                    d_coef_e*((f[idx_f_LLLLL] + f[idx_f_R]    )*(g[idx_g_LLLLL] + g[idx_g_R]    )*(h[idx_h_LLLLL] + h[idx_h_R]    )  +
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
                
                F_face_x[idx_face_x] += dt*quarter*(
                    d_coef_a*((f[idx_f_L]      + f[idx_f_R]     )*(g[idx_g_L]      + g[idx_g_R]     )*(h[idx_h_L]      + h[idx_h_R]     )) +
                    d_coef_b*((f[idx_f_LL]     + f[idx_f_R]     )*(g[idx_g_LL]     + g[idx_g_R]     )*(h[idx_h_LL]     + h[idx_h_R]     )  +
                              (f[idx_f_L]      + f[idx_f_RR]    )*(g[idx_g_L]      + g[idx_g_RR]    )*(h[idx_h_L]      + h[idx_h_RR]    )) +
                    d_coef_c*((f[idx_f_LLL]    + f[idx_f_R]     )*(g[idx_g_LLL]    + g[idx_g_R]     )*(h[idx_h_LLL]    + h[idx_h_R]     )  +
                              (f[idx_f_LL]     + f[idx_f_RR]    )*(g[idx_g_LL]     + g[idx_g_RR]    )*(h[idx_h_LL]     + h[idx_h_RR]    )  +
                              (f[idx_f_L]      + f[idx_f_RRR]   )*(g[idx_g_L]      + g[idx_g_RRR]   )*(h[idx_h_L]      + h[idx_h_RRR]   )) +
                    d_coef_d*((f[idx_f_LLLL]   + f[idx_f_R]     )*(g[idx_g_LLLL]   + g[idx_g_R]     )*(h[idx_h_LLLL]   + h[idx_h_R]     )  +
                              (f[idx_f_LLL]    + f[idx_f_RR]    )*(g[idx_g_LLL]    + g[idx_g_RR]    )*(h[idx_h_LLL]    + h[idx_h_RR]    )  +
                              (f[idx_f_LL]     + f[idx_f_RRR]   )*(g[idx_g_LL]     + g[idx_g_RRR]   )*(h[idx_h_LL]     + h[idx_h_RRR]   )  +
                              (f[idx_f_L]      + f[idx_f_RRRR]  )*(g[idx_g_L]      + g[idx_g_RRRR]  )*(h[idx_h_L]      + h[idx_h_RRRR]  )) +
                    d_coef_e*((f[idx_f_LLLLL]  + f[idx_f_R]     )*(g[idx_g_LLLLL]  + g[idx_g_R]     )*(h[idx_h_LLLLL]  + h[idx_h_R]     )  +
                              (f[idx_f_LLLL]   + f[idx_f_RR]    )*(g[idx_g_LLLL]   + g[idx_g_RR]    )*(h[idx_h_LLLL]   + h[idx_h_RR]    )  +
                              (f[idx_f_LLL]    + f[idx_f_RRR]   )*(g[idx_g_LLL]    + g[idx_g_RRR]   )*(h[idx_h_LLL]    + h[idx_h_RRR]   )  +
                              (f[idx_f_LL]     + f[idx_f_RRRR]  )*(g[idx_g_LL]     + g[idx_g_RRRR]  )*(h[idx_h_LL]     + h[idx_h_RRRR]  )  +
                              (f[idx_f_L]      + f[idx_f_RRRRR] )*(g[idx_g_L]      + g[idx_g_RRRRR] )*(h[idx_h_L]      + h[idx_h_RRRRR] )) +
                    d_coef_f*((f[idx_f_LLLLLL] + f[idx_f_R]     )*(g[idx_g_LLLLLL] + g[idx_g_R]     )*(h[idx_h_LLLLLL] + h[idx_h_R]     )  +
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
                    
                    F_face_x[idx_face_x] += dt*quarter*(
                        d_coef_a*(f[idx_f_L] + f[idx_f_R])*(g[idx_g_L] + g[idx_g_R])*(h[idx_h_L] + h[idx_h_R]));
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
                    
                    F_face_x[idx_face_x] += dt*quarter*(
                        d_coef_a*((f[idx_f_L]  + f[idx_f_R] )*(g[idx_g_L]  + g[idx_g_R] )*(h[idx_h_L]  + h[idx_h_R] )) +
                        d_coef_b*((f[idx_f_LL] + f[idx_f_R] )*(g[idx_g_LL] + g[idx_g_R] )*(h[idx_h_LL] + h[idx_h_R] )  +
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
                    
                    F_face_x[idx_face_x] += dt*quarter*(
                        d_coef_a*((f[idx_f_L]   + f[idx_f_R]  )*(g[idx_g_L]   + g[idx_g_R]  )*(h[idx_h_L]   + h[idx_h_R]  )) +
                        d_coef_b*((f[idx_f_LL]  + f[idx_f_R]  )*(g[idx_g_LL]  + g[idx_g_R]  )*(h[idx_h_LL]  + h[idx_h_R]  )  +
                                  (f[idx_f_L]   + f[idx_f_RR] )*(g[idx_g_L]   + g[idx_g_RR] )*(h[idx_h_L]   + h[idx_h_RR] )) +
                        d_coef_c*((f[idx_f_LLL] + f[idx_f_R]  )*(g[idx_g_LLL] + g[idx_g_R]  )*(h[idx_h_LLL] + h[idx_h_R]  )  +
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
                    
                    F_face_x[idx_face_x] += dt*quarter*(
                        d_coef_a*((f[idx_f_L]    + f[idx_f_R]   )*(g[idx_g_L]    + g[idx_g_R]   )*(h[idx_h_L]    + h[idx_h_R]   )) +
                        d_coef_b*((f[idx_f_LL]   + f[idx_f_R]   )*(g[idx_g_LL]   + g[idx_g_R]   )*(h[idx_h_LL]   + h[idx_h_R]   )  +
                                  (f[idx_f_L]    + f[idx_f_RR]  )*(g[idx_g_L]    + g[idx_g_RR]  )*(h[idx_h_L]    + h[idx_h_RR]  )) +
                        d_coef_c*((f[idx_f_LLL]  + f[idx_f_R]   )*(g[idx_g_LLL]  + g[idx_g_R]   )*(h[idx_h_LLL]  + h[idx_h_R]   )  +
                                  (f[idx_f_LL]   + f[idx_f_RR]  )*(g[idx_g_LL]   + g[idx_g_RR]  )*(h[idx_h_LL]   + h[idx_h_RR]  )  +
                                  (f[idx_f_L]    + f[idx_f_RRR] )*(g[idx_g_L]    + g[idx_g_RRR] )*(h[idx_h_L]    + h[idx_h_RRR] )) +
                        d_coef_d*((f[idx_f_LLLL] + f[idx_f_R]   )*(g[idx_g_LLLL] + g[idx_g_R]   )*(h[idx_h_LLLL] + h[idx_h_R]   )  +
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
                    
                    F_face_x[idx_face_x] += dt*quarter*(
                        d_coef_a*((f[idx_f_L]     + f[idx_f_R]    )*(g[idx_g_L]     + g[idx_g_R]    )*(h[idx_h_L]     + h[idx_h_R]    )) +
                        d_coef_b*((f[idx_f_LL]    + f[idx_f_R]    )*(g[idx_g_LL]    + g[idx_g_R]    )*(h[idx_h_LL]    + h[idx_h_R]    )  +
                                  (f[idx_f_L]     + f[idx_f_RR]   )*(g[idx_g_L]     + g[idx_g_RR]   )*(h[idx_h_L]     + h[idx_h_RR]   )) +
                        d_coef_c*((f[idx_f_LLL]   + f[idx_f_R]    )*(g[idx_g_LLL]   + g[idx_g_R]    )*(h[idx_h_LLL]   + h[idx_h_R]    )  +
                                  (f[idx_f_LL]    + f[idx_f_RR]   )*(g[idx_g_LL]    + g[idx_g_RR]   )*(h[idx_h_LL]    + h[idx_h_RR]   )  +
                                  (f[idx_f_L]     + f[idx_f_RRR]  )*(g[idx_g_L]     + g[idx_g_RRR]  )*(h[idx_h_L]     + h[idx_h_RRR]  )) +
                        d_coef_d*((f[idx_f_LLLL]  + f[idx_f_R]    )*(g[idx_g_LLLL]  + g[idx_g_R]    )*(h[idx_h_LLLL]  + h[idx_h_R]    )  +
                                  (f[idx_f_LLL]   + f[idx_f_RR]   )*(g[idx_g_LLL]   + g[idx_g_RR]   )*(h[idx_h_LLL]   + h[idx_h_RR]   )  +
                                  (f[idx_f_LL]    + f[idx_f_RRR]  )*(g[idx_g_LL]    + g[idx_g_RRR]  )*(h[idx_h_LL]    + h[idx_h_RRR]  )  +
                                  (f[idx_f_L]     + f[idx_f_RRRR] )*(g[idx_g_L]     + g[idx_g_RRRR] )*(h[idx_h_L]     + h[idx_h_RRRR] )) +
                        d_coef_e*((f[idx_f_LLLLL] + f[idx_f_R]    )*(g[idx_g_LLLLL] + g[idx_g_R]    )*(h[idx_h_LLLLL] + h[idx_h_R]    )  +
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
                    
                    F_face_x[idx_face_x] += dt*quarter*(
                        d_coef_a*((f[idx_f_L]      + f[idx_f_R]     )*(g[idx_g_L]      + g[idx_g_R]     )*(h[idx_h_L]      + h[idx_h_R]     )) +
                        d_coef_b*((f[idx_f_LL]     + f[idx_f_R]     )*(g[idx_g_LL]     + g[idx_g_R]     )*(h[idx_h_LL]     + h[idx_h_R]     )  +
                                  (f[idx_f_L]      + f[idx_f_RR]    )*(g[idx_g_L]      + g[idx_g_RR]    )*(h[idx_h_L]      + h[idx_h_RR]    )) +
                        d_coef_c*((f[idx_f_LLL]    + f[idx_f_R]     )*(g[idx_g_LLL]    + g[idx_g_R]     )*(h[idx_h_LLL]    + h[idx_h_R]     )  +
                                  (f[idx_f_LL]     + f[idx_f_RR]    )*(g[idx_g_LL]     + g[idx_g_RR]    )*(h[idx_h_LL]     + h[idx_h_RR]    )  +
                                  (f[idx_f_L]      + f[idx_f_RRR]   )*(g[idx_g_L]      + g[idx_g_RRR]   )*(h[idx_h_L]      + h[idx_h_RRR]   )) +
                        d_coef_d*((f[idx_f_LLLL]   + f[idx_f_R]     )*(g[idx_g_LLLL]   + g[idx_g_R]     )*(h[idx_h_LLLL]   + h[idx_h_R]     )  +
                                  (f[idx_f_LLL]    + f[idx_f_RR]    )*(g[idx_g_LLL]    + g[idx_g_RR]    )*(h[idx_h_LLL]    + h[idx_h_RR]    )  +
                                  (f[idx_f_LL]     + f[idx_f_RRR]   )*(g[idx_g_LL]     + g[idx_g_RRR]   )*(h[idx_h_LL]     + h[idx_h_RRR]   )  +
                                  (f[idx_f_L]      + f[idx_f_RRRR]  )*(g[idx_g_L]      + g[idx_g_RRRR]  )*(h[idx_h_L]      + h[idx_h_RRRR]  )) +
                        d_coef_e*((f[idx_f_LLLLL]  + f[idx_f_R]     )*(g[idx_g_LLLLL]  + g[idx_g_R]     )*(h[idx_h_LLLLL]  + h[idx_h_R]     )  +
                                  (f[idx_f_LLLL]   + f[idx_f_RR]    )*(g[idx_g_LLLL]   + g[idx_g_RR]    )*(h[idx_h_LLLL]   + h[idx_h_RR]    )  +
                                  (f[idx_f_LLL]    + f[idx_f_RRR]   )*(g[idx_g_LLL]    + g[idx_g_RRR]   )*(h[idx_h_LLL]    + h[idx_h_RRR]   )  +
                                  (f[idx_f_LL]     + f[idx_f_RRRR]  )*(g[idx_g_LL]     + g[idx_g_RRRR]  )*(h[idx_h_LL]     + h[idx_h_RRRR]  )  +
                                  (f[idx_f_L]      + f[idx_f_RRRRR] )*(g[idx_g_L]      + g[idx_g_RRRRR] )*(h[idx_h_L]      + h[idx_h_RRRRR] )) +
                        d_coef_f*((f[idx_f_LLLLLL] + f[idx_f_R]     )*(g[idx_g_LLLLLL] + g[idx_g_R]     )*(h[idx_h_LLLLLL] + h[idx_h_R]     )  +
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
                        
                        F_face_x[idx_face_x] += dt*quarter*(
                            d_coef_a*(f[idx_f_L] + f[idx_f_R])*(g[idx_g_L] + g[idx_g_R])*(h[idx_h_L] + h[idx_h_R]));
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
                        
                        F_face_x[idx_face_x] += dt*quarter*(
                            d_coef_a*((f[idx_f_L]  + f[idx_f_R] )*(g[idx_g_L]  + g[idx_g_R] )*(h[idx_h_L]  + h[idx_h_R] )) +
                            d_coef_b*((f[idx_f_LL] + f[idx_f_R] )*(g[idx_g_LL] + g[idx_g_R] )*(h[idx_h_LL] + h[idx_h_R] )  +
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
                        
                        F_face_x[idx_face_x] += dt*quarter*(
                            d_coef_a*((f[idx_f_L]   + f[idx_f_R]  )*(g[idx_g_L]   + g[idx_g_R]  )*(h[idx_h_L]   + h[idx_h_R]  )) +
                            d_coef_b*((f[idx_f_LL]  + f[idx_f_R]  )*(g[idx_g_LL]  + g[idx_g_R]  )*(h[idx_h_LL]  + h[idx_h_R]  )  +
                                      (f[idx_f_L]   + f[idx_f_RR] )*(g[idx_g_L]   + g[idx_g_RR] )*(h[idx_h_L]   + h[idx_h_RR] )) +
                            d_coef_c*((f[idx_f_LLL] + f[idx_f_R]  )*(g[idx_g_LLL] + g[idx_g_R]  )*(h[idx_h_LLL] + h[idx_h_R]  )  +
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
                        
                        F_face_x[idx_face_x] += dt*quarter*(
                            d_coef_a*((f[idx_f_L]    + f[idx_f_R]   )*(g[idx_g_L]    + g[idx_g_R]   )*(h[idx_h_L]    + h[idx_h_R]   )) +
                            d_coef_b*((f[idx_f_LL]   + f[idx_f_R]   )*(g[idx_g_LL]   + g[idx_g_R]   )*(h[idx_h_LL]   + h[idx_h_R]   )  +
                                      (f[idx_f_L]    + f[idx_f_RR]  )*(g[idx_g_L]    + g[idx_g_RR]  )*(h[idx_h_L]    + h[idx_h_RR]  )) +
                            d_coef_c*((f[idx_f_LLL]  + f[idx_f_R]   )*(g[idx_g_LLL]  + g[idx_g_R]   )*(h[idx_h_LLL]  + h[idx_h_R]   )  +
                                      (f[idx_f_LL]   + f[idx_f_RR]  )*(g[idx_g_LL]   + g[idx_g_RR]  )*(h[idx_h_LL]   + h[idx_h_RR]  )  +
                                      (f[idx_f_L]    + f[idx_f_RRR] )*(g[idx_g_L]    + g[idx_g_RRR] )*(h[idx_h_L]    + h[idx_h_RRR] )) +
                            d_coef_d*((f[idx_f_LLLL] + f[idx_f_R]   )*(g[idx_g_LLLL] + g[idx_g_R]   )*(h[idx_h_LLLL] + h[idx_h_R]   )  +
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
                        
                        F_face_x[idx_face_x] += dt*quarter*(
                            d_coef_a*((f[idx_f_L]     + f[idx_f_R]    )*(g[idx_g_L]     + g[idx_g_R]    )*(h[idx_h_L]     + h[idx_h_R]    )) +
                            d_coef_b*((f[idx_f_LL]    + f[idx_f_R]    )*(g[idx_g_LL]    + g[idx_g_R]    )*(h[idx_h_LL]    + h[idx_h_R]    )  +
                                      (f[idx_f_L]     + f[idx_f_RR]   )*(g[idx_g_L]     + g[idx_g_RR]   )*(h[idx_h_L]     + h[idx_h_RR]   )) +
                            d_coef_c*((f[idx_f_LLL]   + f[idx_f_R]    )*(g[idx_g_LLL]   + g[idx_g_R]    )*(h[idx_h_LLL]   + h[idx_h_R]    )  +
                                      (f[idx_f_LL]    + f[idx_f_RR]   )*(g[idx_g_LL]    + g[idx_g_RR]   )*(h[idx_h_LL]    + h[idx_h_RR]   )  +
                                      (f[idx_f_L]     + f[idx_f_RRR]  )*(g[idx_g_L]     + g[idx_g_RRR]  )*(h[idx_h_L]     + h[idx_h_RRR]  )) +
                            d_coef_d*((f[idx_f_LLLL]  + f[idx_f_R]    )*(g[idx_g_LLLL]  + g[idx_g_R]    )*(h[idx_h_LLLL]  + h[idx_h_R]    )  +
                                      (f[idx_f_LLL]   + f[idx_f_RR]   )*(g[idx_g_LLL]   + g[idx_g_RR]   )*(h[idx_h_LLL]   + h[idx_h_RR]   )  +
                                      (f[idx_f_LL]    + f[idx_f_RRR]  )*(g[idx_g_LL]    + g[idx_g_RRR]  )*(h[idx_h_LL]    + h[idx_h_RRR]  )  +
                                      (f[idx_f_L]     + f[idx_f_RRRR] )*(g[idx_g_L]     + g[idx_g_RRRR] )*(h[idx_h_L]     + h[idx_h_RRRR] )) +
                            d_coef_e*((f[idx_f_LLLLL] + f[idx_f_R]    )*(g[idx_g_LLLLL] + g[idx_g_R]    )*(h[idx_h_LLLLL] + h[idx_h_R]    )  +
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
                        
                        F_face_x[idx_face_x] += dt*quarter*(
                            d_coef_a*((f[idx_f_L]      + f[idx_f_R]     )*(g[idx_g_L]      + g[idx_g_R]     )*(h[idx_h_L]      + h[idx_h_R]     )) +
                            d_coef_b*((f[idx_f_LL]     + f[idx_f_R]     )*(g[idx_g_LL]     + g[idx_g_R]     )*(h[idx_h_LL]     + h[idx_h_R]     )  +
                                      (f[idx_f_L]      + f[idx_f_RR]    )*(g[idx_g_L]      + g[idx_g_RR]    )*(h[idx_h_L]      + h[idx_h_RR]    )) +
                            d_coef_c*((f[idx_f_LLL]    + f[idx_f_R]     )*(g[idx_g_LLL]    + g[idx_g_R]     )*(h[idx_h_LLL]    + h[idx_h_R]     )  +
                                      (f[idx_f_LL]     + f[idx_f_RR]    )*(g[idx_g_LL]     + g[idx_g_RR]    )*(h[idx_h_LL]     + h[idx_h_RR]    )  +
                                      (f[idx_f_L]      + f[idx_f_RRR]   )*(g[idx_g_L]      + g[idx_g_RRR]   )*(h[idx_h_L]      + h[idx_h_RRR]   )) +
                            d_coef_d*((f[idx_f_LLLL]   + f[idx_f_R]     )*(g[idx_g_LLLL]   + g[idx_g_R]     )*(h[idx_h_LLLL]   + h[idx_h_R]     )  +
                                      (f[idx_f_LLL]    + f[idx_f_RR]    )*(g[idx_g_LLL]    + g[idx_g_RR]    )*(h[idx_h_LLL]    + h[idx_h_RR]    )  +
                                      (f[idx_f_LL]     + f[idx_f_RRR]   )*(g[idx_g_LL]     + g[idx_g_RRR]   )*(h[idx_h_LL]     + h[idx_h_RRR]   )  +
                                      (f[idx_f_L]      + f[idx_f_RRRR]  )*(g[idx_g_L]      + g[idx_g_RRRR]  )*(h[idx_h_L]      + h[idx_h_RRRR]  )) +
                            d_coef_e*((f[idx_f_LLLLL]  + f[idx_f_R]     )*(g[idx_g_LLLLL]  + g[idx_g_R]     )*(h[idx_h_LLLLL]  + h[idx_h_R]     )  +
                                      (f[idx_f_LLLL]   + f[idx_f_RR]    )*(g[idx_g_LLLL]   + g[idx_g_RR]    )*(h[idx_h_LLLL]   + h[idx_h_RR]    )  +
                                      (f[idx_f_LLL]    + f[idx_f_RRR]   )*(g[idx_g_LLL]    + g[idx_g_RRR]   )*(h[idx_h_LLL]    + h[idx_h_RRR]   )  +
                                      (f[idx_f_LL]     + f[idx_f_RRRR]  )*(g[idx_g_LL]     + g[idx_g_RRRR]  )*(h[idx_h_LL]     + h[idx_h_RRRR]  )  +
                                      (f[idx_f_L]      + f[idx_f_RRRRR] )*(g[idx_g_L]      + g[idx_g_RRRRR] )*(h[idx_h_L]      + h[idx_h_RRRRR] )) +
                            d_coef_f*((f[idx_f_LLLLLL] + f[idx_f_R]     )*(g[idx_g_LLLLLL] + g[idx_g_R]     )*(h[idx_h_LLLLLL] + h[idx_h_R]     )  +
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


/*
 * Add linear term to convective flux in y-direction.
 */
void
ConvectiveFluxReconstructorKEP::addLinearTermToConvectiveFluxY(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_f,
    const int component_idx_flux,
    const int component_idx_f,
    const double dt) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR("ConvectiveFluxReconstructorKEP::addLinearTermToConvectiveFluxY:"
            " Cannot compute linear term in y-direction for 1D problem!");
    }
#endif
    
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
    
    double* F_face_y = data_convective_flux->getPointer(1, component_idx_flux);
    double* f        = data_f->getPointer(component_idx_f);
    
    if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        /*
         * Get the numbers of ghost cells and the dimensions of the ghost cell box.
         */
        
        const int num_ghosts_0_f = num_ghosts_f[0];
        const int num_ghosts_1_f = num_ghosts_f[1];
        const int ghostcell_dim_0_f = ghostcell_dims_f[0];
        
        if (d_stencil_width == 3)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_B = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    F_face_y[idx_face_y] += dt*(
                        d_coef_a*(f[idx_f_B] + f[idx_f_T]));
                }
            }
        }
        else if (d_stencil_width == 5)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BB = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B  = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T  = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    F_face_y[idx_face_y] += dt*(
                        d_coef_a*((f[idx_f_B]  + f[idx_f_T] )) +
                        d_coef_b*((f[idx_f_BB] + f[idx_f_T] )  +
                                  (f[idx_f_B]  + f[idx_f_TT])));
                }
            }
        }
        else if (d_stencil_width == 7)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BBB = (i + num_ghosts_0_f) +
                        (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BB  = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B   = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T   = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT  = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTT = (i + num_ghosts_0_f) +
                        (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    F_face_y[idx_face_y] += dt*(
                        d_coef_a*((f[idx_f_B]   + f[idx_f_T]  )) +
                        d_coef_b*((f[idx_f_BB]  + f[idx_f_T]  )  +
                                  (f[idx_f_B]   + f[idx_f_TT] )) +
                        d_coef_c*((f[idx_f_BBB] + f[idx_f_T]  )  +
                                  (f[idx_f_BB]  + f[idx_f_TT] )  +
                                  (f[idx_f_B]   + f[idx_f_TTT])));
                }
            }
        }
        else if (d_stencil_width == 9)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BBBB = (i + num_ghosts_0_f) +
                        (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBB  = (i + num_ghosts_0_f) +
                        (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BB   = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B    = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T    = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT   = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTT  = (i + num_ghosts_0_f) +
                        (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTT = (i + num_ghosts_0_f) +
                        (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    F_face_y[idx_face_y] += dt*(
                        d_coef_a*((f[idx_f_B]    + f[idx_f_T]   )) +
                        d_coef_b*((f[idx_f_BB]   + f[idx_f_T]   )  +
                                  (f[idx_f_B]    + f[idx_f_TT]  )) +
                        d_coef_c*((f[idx_f_BBB]  + f[idx_f_T]   )  +
                                  (f[idx_f_BB]   + f[idx_f_TT]  )  +
                                  (f[idx_f_B]    + f[idx_f_TTT] )) +
                        d_coef_d*((f[idx_f_BBBB] + f[idx_f_T]   )  +
                                  (f[idx_f_BBB]  + f[idx_f_TT]  )  +
                                  (f[idx_f_BB]   + f[idx_f_TTT] )  +
                                  (f[idx_f_B]    + f[idx_f_TTTT])));
                }
            }
        }
        else if (d_stencil_width == 11)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BBBBB = (i + num_ghosts_0_f) +
                        (j - 5 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBBB  = (i + num_ghosts_0_f) +
                        (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBB   = (i + num_ghosts_0_f) +
                        (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BB    = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B     = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T     = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT    = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTT   = (i + num_ghosts_0_f) +
                        (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTT  = (i + num_ghosts_0_f) +
                        (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTTT = (i + num_ghosts_0_f) +
                        (j + 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    F_face_y[idx_face_y] += dt*(
                        d_coef_a*((f[idx_f_B]     + f[idx_f_T]    )) +
                        d_coef_b*((f[idx_f_BB]    + f[idx_f_T]    )  +
                                  (f[idx_f_B]     + f[idx_f_TT]   )) +
                        d_coef_c*((f[idx_f_BBB]   + f[idx_f_T]    )  +
                                  (f[idx_f_BB]    + f[idx_f_TT]   )  +
                                  (f[idx_f_B]     + f[idx_f_TTT]  )) +
                        d_coef_d*((f[idx_f_BBBB]  + f[idx_f_T]    )  +
                                  (f[idx_f_BBB]   + f[idx_f_TT]   )  +
                                  (f[idx_f_BB]    + f[idx_f_TTT]  )  +
                                  (f[idx_f_B]     + f[idx_f_TTTT] )) +
                        d_coef_e*((f[idx_f_BBBBB] + f[idx_f_T]    )  +
                                  (f[idx_f_BBBB]  + f[idx_f_TT]   )  +
                                  (f[idx_f_BBB]   + f[idx_f_TTT]  )  +
                                  (f[idx_f_BB]    + f[idx_f_TTTT] )  +
                                  (f[idx_f_B]     + f[idx_f_TTTTT])));
                }
            }
        }
        else if (d_stencil_width == 13)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BBBBBB = (i + num_ghosts_0_f) +
                        (j - 6 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBBBB  = (i + num_ghosts_0_f) +
                        (j - 5 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBBB   = (i + num_ghosts_0_f) +
                        (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBB    = (i + num_ghosts_0_f) +
                        (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BB     = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B      = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T      = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT     = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTT    = (i + num_ghosts_0_f) +
                        (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTT   = (i + num_ghosts_0_f) +
                        (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTTT  = (i + num_ghosts_0_f) +
                        (j + 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTTTT = (i + num_ghosts_0_f) +
                        (j + 5 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    F_face_y[idx_face_y] += dt*(
                        d_coef_a*((f[idx_f_B]      + f[idx_f_T]     )) +
                        d_coef_b*((f[idx_f_BB]     + f[idx_f_T]     )  +
                                  (f[idx_f_B]      + f[idx_f_TT]    )) +
                        d_coef_c*((f[idx_f_BBB]    + f[idx_f_T]     )  +
                                  (f[idx_f_BB]     + f[idx_f_TT]    )  +
                                  (f[idx_f_B]      + f[idx_f_TTT]   )) +
                        d_coef_d*((f[idx_f_BBBB]   + f[idx_f_T]     )  +
                                  (f[idx_f_BBB]    + f[idx_f_TT]    )  +
                                  (f[idx_f_BB]     + f[idx_f_TTT]   )  +
                                  (f[idx_f_B]      + f[idx_f_TTTT]  )) +
                        d_coef_e*((f[idx_f_BBBBB]  + f[idx_f_T]     )  +
                                  (f[idx_f_BBBB]   + f[idx_f_TT]    )  +
                                  (f[idx_f_BBB]    + f[idx_f_TTT]   )  +
                                  (f[idx_f_BB]     + f[idx_f_TTTT]  )  +
                                  (f[idx_f_B]      + f[idx_f_TTTTT] )) +
                        d_coef_f*((f[idx_f_BBBBBB] + f[idx_f_T]     )  +
                                  (f[idx_f_BBBBB]  + f[idx_f_TT]    )  +
                                  (f[idx_f_BBBB]   + f[idx_f_TTT]   )  +
                                  (f[idx_f_BBB]    + f[idx_f_TTTT]  )  +
                                  (f[idx_f_BB]     + f[idx_f_TTTTT] )  +
                                  (f[idx_f_B]      + f[idx_f_TTTTTT])));
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
         * Get the numbers of ghost cells and the dimensions of the ghost cell box.
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
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_B = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_y[idx_face_y] += dt*(
                            d_coef_a*(f[idx_f_B] + f[idx_f_T]));
                    }
                }
            }
        }
        else if (d_stencil_width == 5)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BB = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B  = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_y[idx_face_y] += dt*(
                            d_coef_a*((f[idx_f_B]  + f[idx_f_T] )) +
                            d_coef_b*((f[idx_f_BB] + f[idx_f_T] )  +
                                      (f[idx_f_B]  + f[idx_f_TT])));
                    }
                }
            }
        }
        else if (d_stencil_width == 7)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BBB = (i + num_ghosts_0_f) +
                            (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB  = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B   = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT  = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTT = (i + num_ghosts_0_f) +
                            (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_y[idx_face_y] += dt*(
                            d_coef_a*((f[idx_f_B]   + f[idx_f_T]  )) +
                            d_coef_b*((f[idx_f_BB]  + f[idx_f_T]  )  +
                                      (f[idx_f_B]   + f[idx_f_TT] )) +
                            d_coef_c*((f[idx_f_BBB] + f[idx_f_T]  )  +
                                      (f[idx_f_BB]  + f[idx_f_TT] )  +
                                      (f[idx_f_B]   + f[idx_f_TTT])));
                    }
                }
            }
        }
        else if (d_stencil_width == 9)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BBBB = (i + num_ghosts_0_f) +
                            (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB  = (i + num_ghosts_0_f) +
                            (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB   = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B    = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT   = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTT  = (i + num_ghosts_0_f) +
                            (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTT = (i + num_ghosts_0_f) +
                            (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_y[idx_face_y] += dt*(
                            d_coef_a*((f[idx_f_B]    + f[idx_f_T]   )) +
                            d_coef_b*((f[idx_f_BB]   + f[idx_f_T]   )  +
                                      (f[idx_f_B]    + f[idx_f_TT]  )) +
                            d_coef_c*((f[idx_f_BBB]  + f[idx_f_T]   )  +
                                      (f[idx_f_BB]   + f[idx_f_TT]  )  +
                                      (f[idx_f_B]    + f[idx_f_TTT] )) +
                            d_coef_d*((f[idx_f_BBBB] + f[idx_f_T]   )  +
                                      (f[idx_f_BBB]  + f[idx_f_TT]  )  +
                                      (f[idx_f_BB]   + f[idx_f_TTT] )  +
                                      (f[idx_f_B]    + f[idx_f_TTTT])));
                    }
                }
            }
        }
        else if (d_stencil_width == 11)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BBBBB = (i + num_ghosts_0_f) +
                            (j - 5 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBB  = (i + num_ghosts_0_f) +
                            (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB   = (i + num_ghosts_0_f) +
                            (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB    = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B     = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT    = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTT   = (i + num_ghosts_0_f) +
                            (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTT  = (i + num_ghosts_0_f) +
                            (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTTT = (i + num_ghosts_0_f) +
                            (j + 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_y[idx_face_y] += dt*(
                            d_coef_a*((f[idx_f_B]     + f[idx_f_T]    )) +
                            d_coef_b*((f[idx_f_BB]    + f[idx_f_T]    )  +
                                      (f[idx_f_B]     + f[idx_f_TT]   )) +
                            d_coef_c*((f[idx_f_BBB]   + f[idx_f_T]    )  +
                                      (f[idx_f_BB]    + f[idx_f_TT]   )  +
                                      (f[idx_f_B]     + f[idx_f_TTT]  )) +
                            d_coef_d*((f[idx_f_BBBB]  + f[idx_f_T]    )  +
                                      (f[idx_f_BBB]   + f[idx_f_TT]   )  +
                                      (f[idx_f_BB]    + f[idx_f_TTT]  )  +
                                      (f[idx_f_B]     + f[idx_f_TTTT] )) +
                            d_coef_e*((f[idx_f_BBBBB] + f[idx_f_T]    )  +
                                      (f[idx_f_BBBB]  + f[idx_f_TT]   )  +
                                      (f[idx_f_BBB]   + f[idx_f_TTT]  )  +
                                      (f[idx_f_BB]    + f[idx_f_TTTT] )  +
                                      (f[idx_f_B]     + f[idx_f_TTTTT])));
                    }
                }
            }
        }
        else if (d_stencil_width == 13)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BBBBBB = (i + num_ghosts_0_f) +
                            (j - 6 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBBB  = (i + num_ghosts_0_f) +
                            (j - 5 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBB   = (i + num_ghosts_0_f) +
                            (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB    = (i + num_ghosts_0_f) +
                            (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB     = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B      = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T      = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT     = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTT    = (i + num_ghosts_0_f) +
                            (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTT   = (i + num_ghosts_0_f) +
                            (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTTT  = (i + num_ghosts_0_f) +
                            (j + 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTTTT = (i + num_ghosts_0_f) +
                            (j + 5 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_y[idx_face_y] += dt*(
                            d_coef_a*((f[idx_f_B]      + f[idx_f_T]     )) +
                            d_coef_b*((f[idx_f_BB]     + f[idx_f_T]     )  +
                                      (f[idx_f_B]      + f[idx_f_TT]    )) +
                            d_coef_c*((f[idx_f_BBB]    + f[idx_f_T]     )  +
                                      (f[idx_f_BB]     + f[idx_f_TT]    )  +
                                      (f[idx_f_B]      + f[idx_f_TTT]   )) +
                            d_coef_d*((f[idx_f_BBBB]   + f[idx_f_T]     )  +
                                      (f[idx_f_BBB]    + f[idx_f_TT]    )  +
                                      (f[idx_f_BB]     + f[idx_f_TTT]   )  +
                                      (f[idx_f_B]      + f[idx_f_TTTT]  )) +
                            d_coef_e*((f[idx_f_BBBBB]  + f[idx_f_T]     )  +
                                      (f[idx_f_BBBB]   + f[idx_f_TT]    )  +
                                      (f[idx_f_BBB]    + f[idx_f_TTT]   )  +
                                      (f[idx_f_BB]     + f[idx_f_TTTT]  )  +
                                      (f[idx_f_B]      + f[idx_f_TTTTT] )) +
                            d_coef_f*((f[idx_f_BBBBBB] + f[idx_f_T]     )  +
                                      (f[idx_f_BBBBB]  + f[idx_f_TT]    )  +
                                      (f[idx_f_BBBB]   + f[idx_f_TTT]   )  +
                                      (f[idx_f_BBB]    + f[idx_f_TTTT]  )  +
                                      (f[idx_f_BB]     + f[idx_f_TTTTT] )  +
                                      (f[idx_f_B]      + f[idx_f_TTTTTT])));
                    }
                }
            }
        }
    }
}


/*
 * Add quadratic term to convective flux in y-direction.
 */
void
ConvectiveFluxReconstructorKEP::addQuadraticTermToConvectiveFluxY(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_f,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_g,
    const int component_idx_flux,
    const int component_idx_f,
    const int component_idx_g,
    const double dt) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR("ConvectiveFluxReconstructorKEP::addQuadraticTermToConvectiveFluxY:"
            " Cannot compute quadratic term in y-direction for 1D problem!");
    }
#endif
    
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
    
    double* F_face_y = data_convective_flux->getPointer(1, component_idx_flux);
    double* f        = data_f->getPointer(component_idx_f);
    double* g        = data_g->getPointer(component_idx_g);
    
    const double half = double(1)/double(2);
    
    if (d_dim == tbox::Dimension(2))
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
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_B = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_B = (i + num_ghosts_0_g) +
                        (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_T = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    F_face_y[idx_face_y] += dt*half*(
                        d_coef_a*(f[idx_f_B] + f[idx_f_T])*(g[idx_g_B] + g[idx_g_T]));
                }
            }
        }
        else if (d_stencil_width == 5)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BB = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B  = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T  = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_BB = (i + num_ghosts_0_g) +
                        (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_B  = (i + num_ghosts_0_g) +
                        (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_T  = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TT = (i + num_ghosts_0_g) +
                        (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    F_face_y[idx_face_y] += dt*half*(
                        d_coef_a*((f[idx_f_B]  + f[idx_f_T] )*(g[idx_g_B]  + g[idx_g_T] )) +
                        d_coef_b*((f[idx_f_BB] + f[idx_f_T] )*(g[idx_g_BB] + g[idx_g_T] )  +
                                  (f[idx_f_B]  + f[idx_f_TT])*(g[idx_g_B]  + g[idx_g_TT])));
                }
            }
        }
        else if (d_stencil_width == 7)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BBB = (i + num_ghosts_0_f) +
                        (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BB  = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B   = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T   = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT  = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTT = (i + num_ghosts_0_f) +
                        (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_BBB = (i + num_ghosts_0_g) +
                        (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BB  = (i + num_ghosts_0_g) +
                        (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_B   = (i + num_ghosts_0_g) +
                        (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_T   = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TT  = (i + num_ghosts_0_g) +
                        (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTT = (i + num_ghosts_0_g) +
                        (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    F_face_y[idx_face_y] += dt*half*(
                        d_coef_a*((f[idx_f_B]   + f[idx_f_T]  )*(g[idx_g_B]   + g[idx_g_T]  )) +
                        d_coef_b*((f[idx_f_BB]  + f[idx_f_T]  )*(g[idx_g_BB]  + g[idx_g_T]  )  +
                                  (f[idx_f_B]   + f[idx_f_TT] )*(g[idx_g_B]   + g[idx_g_TT] )) +
                        d_coef_c*((f[idx_f_BBB] + f[idx_f_T]  )*(g[idx_g_BBB] + g[idx_g_T]  )  +
                                  (f[idx_f_BB]  + f[idx_f_TT] )*(g[idx_g_BB]  + g[idx_g_TT] )  +
                                  (f[idx_f_B]   + f[idx_f_TTT])*(g[idx_g_B]   + g[idx_g_TTT])));
                }
            }
        }
        else if (d_stencil_width == 9)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BBBB = (i + num_ghosts_0_f) +
                        (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBB  = (i + num_ghosts_0_f) +
                        (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BB   = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B    = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T    = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT   = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTT  = (i + num_ghosts_0_f) +
                        (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTT = (i + num_ghosts_0_f) +
                        (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_BBBB = (i + num_ghosts_0_g) +
                        (j - 4 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BBB  = (i + num_ghosts_0_g) +
                        (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BB   = (i + num_ghosts_0_g) +
                        (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_B    = (i + num_ghosts_0_g) +
                        (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_T    = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TT   = (i + num_ghosts_0_g) +
                        (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTT  = (i + num_ghosts_0_g) +
                        (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTTT = (i + num_ghosts_0_g) +
                        (j + 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    F_face_y[idx_face_y] += dt*half*(
                        d_coef_a*((f[idx_f_B]    + f[idx_f_T]   )*(g[idx_g_B]    + g[idx_g_T]   )) +
                        d_coef_b*((f[idx_f_BB]   + f[idx_f_T]   )*(g[idx_g_BB]   + g[idx_g_T]   )  +
                                  (f[idx_f_B]    + f[idx_f_TT]  )*(g[idx_g_B]    + g[idx_g_TT]  )) +
                        d_coef_c*((f[idx_f_BBB]  + f[idx_f_T]   )*(g[idx_g_BBB]  + g[idx_g_T]   )  +
                                  (f[idx_f_BB]   + f[idx_f_TT]  )*(g[idx_g_BB]   + g[idx_g_TT]  )  +
                                  (f[idx_f_B]    + f[idx_f_TTT] )*(g[idx_g_B]    + g[idx_g_TTT] )) +
                        d_coef_d*((f[idx_f_BBBB] + f[idx_f_T]   )*(g[idx_g_BBBB] + g[idx_g_T]   )  +
                                  (f[idx_f_BBB]  + f[idx_f_TT]  )*(g[idx_g_BBB]  + g[idx_g_TT]  )  +
                                  (f[idx_f_BB]   + f[idx_f_TTT] )*(g[idx_g_BB]   + g[idx_g_TTT] )  +
                                  (f[idx_f_B]    + f[idx_f_TTTT])*(g[idx_g_B]    + g[idx_g_TTTT])));
                }
            }
        }
        else if (d_stencil_width == 11)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BBBBB = (i + num_ghosts_0_f) +
                        (j - 5 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBBB  = (i + num_ghosts_0_f) +
                        (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBB   = (i + num_ghosts_0_f) +
                        (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BB    = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B     = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T     = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT    = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTT   = (i + num_ghosts_0_f) +
                        (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTT  = (i + num_ghosts_0_f) +
                        (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTTT = (i + num_ghosts_0_f) +
                        (j + 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_BBBBB = (i + num_ghosts_0_g) +
                        (j - 5 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BBBB  = (i + num_ghosts_0_g) +
                        (j - 4 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BBB   = (i + num_ghosts_0_g) +
                        (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BB    = (i + num_ghosts_0_g) +
                        (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_B     = (i + num_ghosts_0_g) +
                        (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_T     = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TT    = (i + num_ghosts_0_g) +
                        (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTT   = (i + num_ghosts_0_g) +
                        (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTTT  = (i + num_ghosts_0_g) +
                        (j + 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTTTT = (i + num_ghosts_0_g) +
                        (j + 4 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    F_face_y[idx_face_y] += dt*half*(
                        d_coef_a*((f[idx_f_B]     + f[idx_f_T]    )*(g[idx_g_B]     + g[idx_g_T]    )) +
                        d_coef_b*((f[idx_f_BB]    + f[idx_f_T]    )*(g[idx_g_BB]    + g[idx_g_T]    )  +
                                  (f[idx_f_B]     + f[idx_f_TT]   )*(g[idx_g_B]     + g[idx_g_TT]   )) +
                        d_coef_c*((f[idx_f_BBB]   + f[idx_f_T]    )*(g[idx_g_BBB]   + g[idx_g_T]    )  +
                                  (f[idx_f_BB]    + f[idx_f_TT]   )*(g[idx_g_BB]    + g[idx_g_TT]   )  +
                                  (f[idx_f_B]     + f[idx_f_TTT]  )*(g[idx_g_B]     + g[idx_g_TTT]  )) +
                        d_coef_d*((f[idx_f_BBBB]  + f[idx_f_T]    )*(g[idx_g_BBBB]  + g[idx_g_T]    )  +
                                  (f[idx_f_BBB]   + f[idx_f_TT]   )*(g[idx_g_BBB]   + g[idx_g_TT]   )  +
                                  (f[idx_f_BB]    + f[idx_f_TTT]  )*(g[idx_g_BB]    + g[idx_g_TTT]  )  +
                                  (f[idx_f_B]     + f[idx_f_TTTT] )*(g[idx_g_B]     + g[idx_g_TTTT] )) +
                        d_coef_e*((f[idx_f_BBBBB] + f[idx_f_T]    )*(g[idx_g_BBBBB] + g[idx_g_T]    )  +
                                  (f[idx_f_BBBB]  + f[idx_f_TT]   )*(g[idx_g_BBBB]  + g[idx_g_TT]   )  +
                                  (f[idx_f_BBB]   + f[idx_f_TTT]  )*(g[idx_g_BBB]   + g[idx_g_TTT]  )  +
                                  (f[idx_f_BB]    + f[idx_f_TTTT] )*(g[idx_g_BB]    + g[idx_g_TTTT] )  +
                                  (f[idx_f_B]     + f[idx_f_TTTTT])*(g[idx_g_B]     + g[idx_g_TTTTT])));
                }
            }
        }
        else if (d_stencil_width == 13)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BBBBBB = (i + num_ghosts_0_f) +
                        (j - 6 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBBBB  = (i + num_ghosts_0_f) +
                        (j - 5 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBBB   = (i + num_ghosts_0_f) +
                        (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBB    = (i + num_ghosts_0_f) +
                        (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BB     = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B      = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T      = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT     = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTT    = (i + num_ghosts_0_f) +
                        (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTT   = (i + num_ghosts_0_f) +
                        (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTTT  = (i + num_ghosts_0_f) +
                        (j + 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTTTT = (i + num_ghosts_0_f) +
                        (j + 5 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_BBBBBB = (i + num_ghosts_0_g) +
                        (j - 6 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BBBBB  = (i + num_ghosts_0_g) +
                        (j - 5 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BBBB   = (i + num_ghosts_0_g) +
                        (j - 4 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BBB    = (i + num_ghosts_0_g) +
                        (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BB     = (i + num_ghosts_0_g) +
                        (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_B      = (i + num_ghosts_0_g) +
                        (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_T      = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TT     = (i + num_ghosts_0_g) +
                        (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTT    = (i + num_ghosts_0_g) +
                        (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTTT   = (i + num_ghosts_0_g) +
                        (j + 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTTTT  = (i + num_ghosts_0_g) +
                        (j + 4 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTTTTT = (i + num_ghosts_0_g) +
                        (j + 5 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    F_face_y[idx_face_y] += dt*half*(
                        d_coef_a*((f[idx_f_B]      + f[idx_f_T]     )*(g[idx_g_B]      + g[idx_g_T]     )) +
                        d_coef_b*((f[idx_f_BB]     + f[idx_f_T]     )*(g[idx_g_BB]     + g[idx_g_T]     )  +
                                  (f[idx_f_B]      + f[idx_f_TT]    )*(g[idx_g_B]      + g[idx_g_TT]    )) +
                        d_coef_c*((f[idx_f_BBB]    + f[idx_f_T]     )*(g[idx_g_BBB]    + g[idx_g_T]     )  +
                                  (f[idx_f_BB]     + f[idx_f_TT]    )*(g[idx_g_BB]     + g[idx_g_TT]    )  +
                                  (f[idx_f_B]      + f[idx_f_TTT]   )*(g[idx_g_B]      + g[idx_g_TTT]   )) +
                        d_coef_d*((f[idx_f_BBBB]   + f[idx_f_T]     )*(g[idx_g_BBBB]   + g[idx_g_T]     )  +
                                  (f[idx_f_BBB]    + f[idx_f_TT]    )*(g[idx_g_BBB]    + g[idx_g_TT]    )  +
                                  (f[idx_f_BB]     + f[idx_f_TTT]   )*(g[idx_g_BB]     + g[idx_g_TTT]   )  +
                                  (f[idx_f_B]      + f[idx_f_TTTT]  )*(g[idx_g_B]      + g[idx_g_TTTT]  )) +
                        d_coef_e*((f[idx_f_BBBBB]  + f[idx_f_T]     )*(g[idx_g_BBBBB]  + g[idx_g_T]     )  +
                                  (f[idx_f_BBBB]   + f[idx_f_TT]    )*(g[idx_g_BBBB]   + g[idx_g_TT]    )  +
                                  (f[idx_f_BBB]    + f[idx_f_TTT]   )*(g[idx_g_BBB]    + g[idx_g_TTT]   )  +
                                  (f[idx_f_BB]     + f[idx_f_TTTT]  )*(g[idx_g_BB]     + g[idx_g_TTTT]  )  +
                                  (f[idx_f_B]      + f[idx_f_TTTTT] )*(g[idx_g_B]      + g[idx_g_TTTTT] )) +
                        d_coef_f*((f[idx_f_BBBBBB] + f[idx_f_T]     )*(g[idx_g_BBBBBB] + g[idx_g_T]     )  +
                                  (f[idx_f_BBBBB]  + f[idx_f_TT]    )*(g[idx_g_BBBBB]  + g[idx_g_TT]    )  +
                                  (f[idx_f_BBBB]   + f[idx_f_TTT]   )*(g[idx_g_BBBB]   + g[idx_g_TTT]   )  +
                                  (f[idx_f_BBB]    + f[idx_f_TTTT]  )*(g[idx_g_BBB]    + g[idx_g_TTTT]  )  +
                                  (f[idx_f_BB]     + f[idx_f_TTTTT] )*(g[idx_g_BB]     + g[idx_g_TTTTT] )  +
                                  (f[idx_f_B]      + f[idx_f_TTTTTT])*(g[idx_g_B]      + g[idx_g_TTTTTT])));
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
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_B = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_B = (i + num_ghosts_0_g) +
                            (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_T = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_y[idx_face_y] += dt*half*(
                            d_coef_a*(f[idx_f_B] + f[idx_f_T])*(g[idx_g_B] + g[idx_g_T]));
                    }
                }
            }
        }
        else if (d_stencil_width == 5)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BB = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B  = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BB = (i + num_ghosts_0_g) +
                            (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B  = (i + num_ghosts_0_g) +
                            (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_T  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TT = (i + num_ghosts_0_g) +
                            (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_y[idx_face_y] += dt*half*(
                            d_coef_a*((f[idx_f_B]  + f[idx_f_T] )*(g[idx_g_B]  + g[idx_g_T] )) +
                            d_coef_b*((f[idx_f_BB] + f[idx_f_T] )*(g[idx_g_BB] + g[idx_g_T] )  +
                                      (f[idx_f_B]  + f[idx_f_TT])*(g[idx_g_B]  + g[idx_g_TT])));
                    }
                }
            }
        }
        else if (d_stencil_width == 7)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BBB = (i + num_ghosts_0_f) +
                            (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB  = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B   = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT  = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTT = (i + num_ghosts_0_f) +
                            (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBB = (i + num_ghosts_0_g) +
                            (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB  = (i + num_ghosts_0_g) +
                            (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B   = (i + num_ghosts_0_g) +
                            (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_T   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TT  = (i + num_ghosts_0_g) +
                            (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTT = (i + num_ghosts_0_g) +
                            (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_y[idx_face_y] += dt*half*(
                            d_coef_a*((f[idx_f_B]   + f[idx_f_T]  )*(g[idx_g_B]   + g[idx_g_T]  )) +
                            d_coef_b*((f[idx_f_BB]  + f[idx_f_T]  )*(g[idx_g_BB]  + g[idx_g_T]  )  +
                                      (f[idx_f_B]   + f[idx_f_TT] )*(g[idx_g_B]   + g[idx_g_TT] )) +
                            d_coef_c*((f[idx_f_BBB] + f[idx_f_T]  )*(g[idx_g_BBB] + g[idx_g_T]  )  +
                                      (f[idx_f_BB]  + f[idx_f_TT] )*(g[idx_g_BB]  + g[idx_g_TT] )  +
                                      (f[idx_f_B]   + f[idx_f_TTT])*(g[idx_g_B]   + g[idx_g_TTT])));
                    }
                }
            }
        }
        else if (d_stencil_width == 9)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BBBB = (i + num_ghosts_0_f) +
                            (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB  = (i + num_ghosts_0_f) +
                            (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB   = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B    = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT   = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTT  = (i + num_ghosts_0_f) +
                            (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTT = (i + num_ghosts_0_f) +
                            (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBBB = (i + num_ghosts_0_g) +
                            (j - 4 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBB  = (i + num_ghosts_0_g) +
                            (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB   = (i + num_ghosts_0_g) +
                            (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B    = (i + num_ghosts_0_g) +
                            (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_T    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TT   = (i + num_ghosts_0_g) +
                            (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTT  = (i + num_ghosts_0_g) +
                            (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTTT = (i + num_ghosts_0_g) +
                            (j + 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_y[idx_face_y] += dt*half*(
                            d_coef_a*((f[idx_f_B]    + f[idx_f_T]   )*(g[idx_g_B]    + g[idx_g_T]   )) +
                            d_coef_b*((f[idx_f_BB]   + f[idx_f_T]   )*(g[idx_g_BB]   + g[idx_g_T]   )  +
                                      (f[idx_f_B]    + f[idx_f_TT]  )*(g[idx_g_B]    + g[idx_g_TT]  )) +
                            d_coef_c*((f[idx_f_BBB]  + f[idx_f_T]   )*(g[idx_g_BBB]  + g[idx_g_T]   )  +
                                      (f[idx_f_BB]   + f[idx_f_TT]  )*(g[idx_g_BB]   + g[idx_g_TT]  )  +
                                      (f[idx_f_B]    + f[idx_f_TTT] )*(g[idx_g_B]    + g[idx_g_TTT] )) +
                            d_coef_d*((f[idx_f_BBBB] + f[idx_f_T]   )*(g[idx_g_BBBB] + g[idx_g_T]   )  +
                                      (f[idx_f_BBB]  + f[idx_f_TT]  )*(g[idx_g_BBB]  + g[idx_g_TT]  )  +
                                      (f[idx_f_BB]   + f[idx_f_TTT] )*(g[idx_g_BB]   + g[idx_g_TTT] )  +
                                      (f[idx_f_B]    + f[idx_f_TTTT])*(g[idx_g_B]    + g[idx_g_TTTT])));
                    }
                }
            }
        }
        else if (d_stencil_width == 11)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BBBBB = (i + num_ghosts_0_f) +
                            (j - 5 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBB  = (i + num_ghosts_0_f) +
                            (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB   = (i + num_ghosts_0_f) +
                            (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB    = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B     = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT    = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTT   = (i + num_ghosts_0_f) +
                            (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTT  = (i + num_ghosts_0_f) +
                            (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTTT = (i + num_ghosts_0_f) +
                            (j + 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBBBB = (i + num_ghosts_0_g) +
                            (j - 5 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBBB  = (i + num_ghosts_0_g) +
                            (j - 4 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBB   = (i + num_ghosts_0_g) +
                            (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB    = (i + num_ghosts_0_g) +
                            (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B     = (i + num_ghosts_0_g) +
                            (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_T     = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TT    = (i + num_ghosts_0_g) +
                            (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTT   = (i + num_ghosts_0_g) +
                            (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTTT  = (i + num_ghosts_0_g) +
                            (j + 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTTTT = (i + num_ghosts_0_g) +
                            (j + 4 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_y[idx_face_y] += dt*half*(
                            d_coef_a*((f[idx_f_B]     + f[idx_f_T]    )*(g[idx_g_B]     + g[idx_g_T]    )) +
                            d_coef_b*((f[idx_f_BB]    + f[idx_f_T]    )*(g[idx_g_BB]    + g[idx_g_T]    )  +
                                      (f[idx_f_B]     + f[idx_f_TT]   )*(g[idx_g_B]     + g[idx_g_TT]   )) +
                            d_coef_c*((f[idx_f_BBB]   + f[idx_f_T]    )*(g[idx_g_BBB]   + g[idx_g_T]    )  +
                                      (f[idx_f_BB]    + f[idx_f_TT]   )*(g[idx_g_BB]    + g[idx_g_TT]   )  +
                                      (f[idx_f_B]     + f[idx_f_TTT]  )*(g[idx_g_B]     + g[idx_g_TTT]  )) +
                            d_coef_d*((f[idx_f_BBBB]  + f[idx_f_T]    )*(g[idx_g_BBBB]  + g[idx_g_T]    )  +
                                      (f[idx_f_BBB]   + f[idx_f_TT]   )*(g[idx_g_BBB]   + g[idx_g_TT]   )  +
                                      (f[idx_f_BB]    + f[idx_f_TTT]  )*(g[idx_g_BB]    + g[idx_g_TTT]  )  +
                                      (f[idx_f_B]     + f[idx_f_TTTT] )*(g[idx_g_B]     + g[idx_g_TTTT] )) +
                            d_coef_e*((f[idx_f_BBBBB] + f[idx_f_T]    )*(g[idx_g_BBBBB] + g[idx_g_T]    )  +
                                      (f[idx_f_BBBB]  + f[idx_f_TT]   )*(g[idx_g_BBBB]  + g[idx_g_TT]   )  +
                                      (f[idx_f_BBB]   + f[idx_f_TTT]  )*(g[idx_g_BBB]   + g[idx_g_TTT]  )  +
                                      (f[idx_f_BB]    + f[idx_f_TTTT] )*(g[idx_g_BB]    + g[idx_g_TTTT] )  +
                                      (f[idx_f_B]     + f[idx_f_TTTTT])*(g[idx_g_B]     + g[idx_g_TTTTT])));
                    }
                }
            }
        }
        else if (d_stencil_width == 13)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BBBBBB = (i + num_ghosts_0_f) +
                            (j - 6 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBBB  = (i + num_ghosts_0_f) +
                            (j - 5 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBB   = (i + num_ghosts_0_f) +
                            (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB    = (i + num_ghosts_0_f) +
                            (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB     = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B      = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T      = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT     = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTT    = (i + num_ghosts_0_f) +
                            (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTT   = (i + num_ghosts_0_f) +
                            (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTTT  = (i + num_ghosts_0_f) +
                            (j + 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTTTT = (i + num_ghosts_0_f) +
                            (j + 5 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBBBBB = (i + num_ghosts_0_g) +
                            (j - 6 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBBBB  = (i + num_ghosts_0_g) +
                            (j - 5 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBBB   = (i + num_ghosts_0_g) +
                            (j - 4 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBB    = (i + num_ghosts_0_g) +
                            (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB     = (i + num_ghosts_0_g) +
                            (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B      = (i + num_ghosts_0_g) +
                            (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_T      = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TT     = (i + num_ghosts_0_g) +
                            (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTT    = (i + num_ghosts_0_g) +
                            (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTTT   = (i + num_ghosts_0_g) +
                            (j + 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTTTT  = (i + num_ghosts_0_g) +
                            (j + 4 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTTTTT = (i + num_ghosts_0_g) +
                            (j + 5 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_y[idx_face_y] += dt*half*(
                            d_coef_a*((f[idx_f_B]      + f[idx_f_T]     )*(g[idx_g_B]      + g[idx_g_T]     )) +
                            d_coef_b*((f[idx_f_BB]     + f[idx_f_T]     )*(g[idx_g_BB]     + g[idx_g_T]     )  +
                                      (f[idx_f_B]      + f[idx_f_TT]    )*(g[idx_g_B]      + g[idx_g_TT]    )) +
                            d_coef_c*((f[idx_f_BBB]    + f[idx_f_T]     )*(g[idx_g_BBB]    + g[idx_g_T]     )  +
                                      (f[idx_f_BB]     + f[idx_f_TT]    )*(g[idx_g_BB]     + g[idx_g_TT]    )  +
                                      (f[idx_f_B]      + f[idx_f_TTT]   )*(g[idx_g_B]      + g[idx_g_TTT]   )) +
                            d_coef_d*((f[idx_f_BBBB]   + f[idx_f_T]     )*(g[idx_g_BBBB]   + g[idx_g_T]     )  +
                                      (f[idx_f_BBB]    + f[idx_f_TT]    )*(g[idx_g_BBB]    + g[idx_g_TT]    )  +
                                      (f[idx_f_BB]     + f[idx_f_TTT]   )*(g[idx_g_BB]     + g[idx_g_TTT]   )  +
                                      (f[idx_f_B]      + f[idx_f_TTTT]  )*(g[idx_g_B]      + g[idx_g_TTTT]  )) +
                            d_coef_e*((f[idx_f_BBBBB]  + f[idx_f_T]     )*(g[idx_g_BBBBB]  + g[idx_g_T]     )  +
                                      (f[idx_f_BBBB]   + f[idx_f_TT]    )*(g[idx_g_BBBB]   + g[idx_g_TT]    )  +
                                      (f[idx_f_BBB]    + f[idx_f_TTT]   )*(g[idx_g_BBB]    + g[idx_g_TTT]   )  +
                                      (f[idx_f_BB]     + f[idx_f_TTTT]  )*(g[idx_g_BB]     + g[idx_g_TTTT]  )  +
                                      (f[idx_f_B]      + f[idx_f_TTTTT] )*(g[idx_g_B]      + g[idx_g_TTTTT] )) +
                            d_coef_f*((f[idx_f_BBBBBB] + f[idx_f_T]     )*(g[idx_g_BBBBBB] + g[idx_g_T]     )  +
                                      (f[idx_f_BBBBB]  + f[idx_f_TT]    )*(g[idx_g_BBBBB]  + g[idx_g_TT]    )  +
                                      (f[idx_f_BBBB]   + f[idx_f_TTT]   )*(g[idx_g_BBBB]   + g[idx_g_TTT]   )  +
                                      (f[idx_f_BBB]    + f[idx_f_TTTT]  )*(g[idx_g_BBB]    + g[idx_g_TTTT]  )  +
                                      (f[idx_f_BB]     + f[idx_f_TTTTT] )*(g[idx_g_BB]     + g[idx_g_TTTTT] )  +
                                      (f[idx_f_B]      + f[idx_f_TTTTTT])*(g[idx_g_B]      + g[idx_g_TTTTTT])));
                    }
                }
            }
        }
    }
}


/*
 * Add cubic term to convective flux in y-direction.
 */
void
ConvectiveFluxReconstructorKEP::addCubicTermToConvectiveFluxY(
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
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR("ConvectiveFluxReconstructorKEP::addCubicTermToConvectiveFluxY:"
            " Cannot compute cubic term in y-direction for 1D problem!");
    }
#endif
    
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
    
    double* F_face_y = data_convective_flux->getPointer(1, component_idx_flux);
    double* f        = data_f->getPointer(component_idx_f);
    double* g        = data_g->getPointer(component_idx_g);
    double* h        = data_h->getPointer(component_idx_h);
    
    const double quarter = double(1)/double(4);
    
    if (d_dim == tbox::Dimension(2))
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
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_B = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_B = (i + num_ghosts_0_g) +
                        (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_T = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_h_B = (i + num_ghosts_0_h) +
                        (j - 1 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_T = (i + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    F_face_y[idx_face_y] += dt*quarter*(
                        d_coef_a*(f[idx_f_B] + f[idx_f_T])*(g[idx_g_B] + g[idx_g_T])*(h[idx_h_B] + h[idx_h_T]));
                }
            }
        }
        else if (d_stencil_width == 5)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BB = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B  = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T  = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_BB = (i + num_ghosts_0_g) +
                        (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_B  = (i + num_ghosts_0_g) +
                        (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_T  = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TT = (i + num_ghosts_0_g) +
                        (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_h_BB = (i + num_ghosts_0_h) +
                        (j - 2 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_B  = (i + num_ghosts_0_h) +
                        (j - 1 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_T  = (i + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TT = (i + num_ghosts_0_h) +
                        (j + 1 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    F_face_y[idx_face_y] += dt*quarter*(
                        d_coef_a*((f[idx_f_B]  + f[idx_f_T] )*(g[idx_g_B]  + g[idx_g_T] )*(h[idx_h_B]  + h[idx_h_T] )) +
                        d_coef_b*((f[idx_f_BB] + f[idx_f_T] )*(g[idx_g_BB] + g[idx_g_T] )*(h[idx_h_BB] + h[idx_h_T] )  +
                                  (f[idx_f_B]  + f[idx_f_TT])*(g[idx_g_B]  + g[idx_g_TT])*(h[idx_h_B]  + h[idx_h_TT])));
                }
            }
        }
        else if (d_stencil_width == 7)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BBB = (i + num_ghosts_0_f) +
                        (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BB  = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B   = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T   = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT  = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTT = (i + num_ghosts_0_f) +
                        (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_BBB = (i + num_ghosts_0_g) +
                        (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BB  = (i + num_ghosts_0_g) +
                        (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_B   = (i + num_ghosts_0_g) +
                        (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_T   = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TT  = (i + num_ghosts_0_g) +
                        (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTT = (i + num_ghosts_0_g) +
                        (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_h_BBB = (i + num_ghosts_0_h) +
                        (j - 3 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_BB  = (i + num_ghosts_0_h) +
                        (j - 2 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_B   = (i + num_ghosts_0_h) +
                        (j - 1 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_T   = (i + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TT  = (i + num_ghosts_0_h) +
                        (j + 1 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TTT = (i + num_ghosts_0_h) +
                        (j + 2 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    F_face_y[idx_face_y] += dt*quarter*(
                        d_coef_a*((f[idx_f_B]   + f[idx_f_T]  )*(g[idx_g_B]   + g[idx_g_T]  )*(h[idx_h_B]   + h[idx_h_T]  )) +
                        d_coef_b*((f[idx_f_BB]  + f[idx_f_T]  )*(g[idx_g_BB]  + g[idx_g_T]  )*(h[idx_h_BB]  + h[idx_h_T]  )  +
                                  (f[idx_f_B]   + f[idx_f_TT] )*(g[idx_g_B]   + g[idx_g_TT] )*(h[idx_h_B]   + h[idx_h_TT] )) +
                        d_coef_c*((f[idx_f_BBB] + f[idx_f_T]  )*(g[idx_g_BBB] + g[idx_g_T]  )*(h[idx_h_BBB] + h[idx_h_T]  )  +
                                  (f[idx_f_BB]  + f[idx_f_TT] )*(g[idx_g_BB]  + g[idx_g_TT] )*(h[idx_h_BB]  + h[idx_h_TT] )  +
                                  (f[idx_f_B]   + f[idx_f_TTT])*(g[idx_g_B]   + g[idx_g_TTT])*(h[idx_h_B]   + h[idx_h_TTT])));
                }
            }
        }
        else if (d_stencil_width == 9)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BBBB = (i + num_ghosts_0_f) +
                        (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBB  = (i + num_ghosts_0_f) +
                        (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BB   = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B    = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T    = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT   = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTT  = (i + num_ghosts_0_f) +
                        (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTT = (i + num_ghosts_0_f) +
                        (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_BBBB = (i + num_ghosts_0_g) +
                        (j - 4 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BBB  = (i + num_ghosts_0_g) +
                        (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BB   = (i + num_ghosts_0_g) +
                        (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_B    = (i + num_ghosts_0_g) +
                        (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_T    = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TT   = (i + num_ghosts_0_g) +
                        (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTT  = (i + num_ghosts_0_g) +
                        (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTTT = (i + num_ghosts_0_g) +
                        (j + 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_h_BBBB = (i + num_ghosts_0_h) +
                        (j - 4 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_BBB  = (i + num_ghosts_0_h) +
                        (j - 3 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_BB   = (i + num_ghosts_0_h) +
                        (j - 2 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_B    = (i + num_ghosts_0_h) +
                        (j - 1 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_T    = (i + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TT   = (i + num_ghosts_0_h) +
                        (j + 1 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TTT  = (i + num_ghosts_0_h) +
                        (j + 2 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TTTT = (i + num_ghosts_0_h) +
                        (j + 3 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    F_face_y[idx_face_y] += dt*quarter*(
                        d_coef_a*((f[idx_f_B]    + f[idx_f_T]   )*(g[idx_g_B]    + g[idx_g_T]   )*(h[idx_h_B]    + h[idx_h_T]   )) +
                        d_coef_b*((f[idx_f_BB]   + f[idx_f_T]   )*(g[idx_g_BB]   + g[idx_g_T]   )*(h[idx_h_BB]   + h[idx_h_T]   )  +
                                  (f[idx_f_B]    + f[idx_f_TT]  )*(g[idx_g_B]    + g[idx_g_TT]  )*(h[idx_h_B]    + h[idx_h_TT]  )) +
                        d_coef_c*((f[idx_f_BBB]  + f[idx_f_T]   )*(g[idx_g_BBB]  + g[idx_g_T]   )*(h[idx_h_BBB]  + h[idx_h_T]   )  +
                                  (f[idx_f_BB]   + f[idx_f_TT]  )*(g[idx_g_BB]   + g[idx_g_TT]  )*(h[idx_h_BB]   + h[idx_h_TT]  )  +
                                  (f[idx_f_B]    + f[idx_f_TTT] )*(g[idx_g_B]    + g[idx_g_TTT] )*(h[idx_h_B]    + h[idx_h_TTT] )) +
                        d_coef_d*((f[idx_f_BBBB] + f[idx_f_T]   )*(g[idx_g_BBBB] + g[idx_g_T]   )*(h[idx_h_BBBB] + h[idx_h_T]   )  +
                                  (f[idx_f_BBB]  + f[idx_f_TT]  )*(g[idx_g_BBB]  + g[idx_g_TT]  )*(h[idx_h_BBB]  + h[idx_h_TT]  )  +
                                  (f[idx_f_BB]   + f[idx_f_TTT] )*(g[idx_g_BB]   + g[idx_g_TTT] )*(h[idx_h_BB]   + h[idx_h_TTT] )  +
                                  (f[idx_f_B]    + f[idx_f_TTTT])*(g[idx_g_B]    + g[idx_g_TTTT])*(h[idx_h_B]    + h[idx_h_TTTT])));
                }
            }
        }
        else if (d_stencil_width == 11)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BBBBB = (i + num_ghosts_0_f) +
                        (j - 5 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBBB  = (i + num_ghosts_0_f) +
                        (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBB   = (i + num_ghosts_0_f) +
                        (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BB    = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B     = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T     = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT    = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTT   = (i + num_ghosts_0_f) +
                        (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTT  = (i + num_ghosts_0_f) +
                        (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTTT = (i + num_ghosts_0_f) +
                        (j + 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_BBBBB = (i + num_ghosts_0_g) +
                        (j - 5 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BBBB  = (i + num_ghosts_0_g) +
                        (j - 4 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BBB   = (i + num_ghosts_0_g) +
                        (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BB    = (i + num_ghosts_0_g) +
                        (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_B     = (i + num_ghosts_0_g) +
                        (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_T     = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TT    = (i + num_ghosts_0_g) +
                        (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTT   = (i + num_ghosts_0_g) +
                        (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTTT  = (i + num_ghosts_0_g) +
                        (j + 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTTTT = (i + num_ghosts_0_g) +
                        (j + 4 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_h_BBBBB = (i + num_ghosts_0_h) +
                        (j - 5 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_BBBB  = (i + num_ghosts_0_h) +
                        (j - 4 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_BBB   = (i + num_ghosts_0_h) +
                        (j - 3 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_BB    = (i + num_ghosts_0_h) +
                        (j - 2 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_B     = (i + num_ghosts_0_h) +
                        (j - 1 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_T     = (i + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TT    = (i + num_ghosts_0_h) +
                        (j + 1 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TTT   = (i + num_ghosts_0_h) +
                        (j + 2 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TTTT  = (i + num_ghosts_0_h) +
                        (j + 3 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TTTTT = (i + num_ghosts_0_h) +
                        (j + 4 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    F_face_y[idx_face_y] += dt*quarter*(
                        d_coef_a*((f[idx_f_B]     + f[idx_f_T]    )*(g[idx_g_B]     + g[idx_g_T]    )*(h[idx_h_B]     + h[idx_h_T]    )) +
                        d_coef_b*((f[idx_f_BB]    + f[idx_f_T]    )*(g[idx_g_BB]    + g[idx_g_T]    )*(h[idx_h_BB]    + h[idx_h_T]    )  +
                                  (f[idx_f_B]     + f[idx_f_TT]   )*(g[idx_g_B]     + g[idx_g_TT]   )*(h[idx_h_B]     + h[idx_h_TT]   )) +
                        d_coef_c*((f[idx_f_BBB]   + f[idx_f_T]    )*(g[idx_g_BBB]   + g[idx_g_T]    )*(h[idx_h_BBB]   + h[idx_h_T]    )  +
                                  (f[idx_f_BB]    + f[idx_f_TT]   )*(g[idx_g_BB]    + g[idx_g_TT]   )*(h[idx_h_BB]    + h[idx_h_TT]   )  +
                                  (f[idx_f_B]     + f[idx_f_TTT]  )*(g[idx_g_B]     + g[idx_g_TTT]  )*(h[idx_h_B]     + h[idx_h_TTT]  )) +
                        d_coef_d*((f[idx_f_BBBB]  + f[idx_f_T]    )*(g[idx_g_BBBB]  + g[idx_g_T]    )*(h[idx_h_BBBB]  + h[idx_h_T]    )  +
                                  (f[idx_f_BBB]   + f[idx_f_TT]   )*(g[idx_g_BBB]   + g[idx_g_TT]   )*(h[idx_h_BBB]   + h[idx_h_TT]   )  +
                                  (f[idx_f_BB]    + f[idx_f_TTT]  )*(g[idx_g_BB]    + g[idx_g_TTT]  )*(h[idx_h_BB]    + h[idx_h_TTT]  )  +
                                  (f[idx_f_B]     + f[idx_f_TTTT] )*(g[idx_g_B]     + g[idx_g_TTTT] )*(h[idx_h_B]     + h[idx_h_TTTT] )) +
                        d_coef_e*((f[idx_f_BBBBB] + f[idx_f_T]    )*(g[idx_g_BBBBB] + g[idx_g_T]    )*(h[idx_h_BBBBB] + h[idx_h_T]    )  +
                                  (f[idx_f_BBBB]  + f[idx_f_TT]   )*(g[idx_g_BBBB]  + g[idx_g_TT]   )*(h[idx_h_BBBB]  + h[idx_h_TT]   )  +
                                  (f[idx_f_BBB]   + f[idx_f_TTT]  )*(g[idx_g_BBB]   + g[idx_g_TTT]  )*(h[idx_h_BBB]   + h[idx_h_TTT]  )  +
                                  (f[idx_f_BB]    + f[idx_f_TTTT] )*(g[idx_g_BB]    + g[idx_g_TTTT] )*(h[idx_h_BB]    + h[idx_h_TTTT] )  +
                                  (f[idx_f_B]     + f[idx_f_TTTTT])*(g[idx_g_B]     + g[idx_g_TTTTT])*(h[idx_h_B]     + h[idx_h_TTTTT])));
                }
            }
        }
        else if (d_stencil_width == 13)
        {
            for (int j = 0; j < interior_dim_1 + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    const int idx_face_y = i +
                        j*interior_dim_0;
                    
                    const int idx_f_BBBBBB = (i + num_ghosts_0_f) +
                        (j - 6 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBBBB  = (i + num_ghosts_0_f) +
                        (j - 5 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBBB   = (i + num_ghosts_0_f) +
                        (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BBB    = (i + num_ghosts_0_f) +
                        (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_BB     = (i + num_ghosts_0_f) +
                        (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_B      = (i + num_ghosts_0_f) +
                        (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_T      = (i + num_ghosts_0_f) +
                        (j + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TT     = (i + num_ghosts_0_f) +
                        (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTT    = (i + num_ghosts_0_f) +
                        (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTT   = (i + num_ghosts_0_f) +
                        (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTTT  = (i + num_ghosts_0_f) +
                        (j + 4 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_f_TTTTTT = (i + num_ghosts_0_f) +
                        (j + 5 + num_ghosts_1_f)*ghostcell_dim_0_f;
                    
                    const int idx_g_BBBBBB = (i + num_ghosts_0_g) +
                        (j - 6 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BBBBB  = (i + num_ghosts_0_g) +
                        (j - 5 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BBBB   = (i + num_ghosts_0_g) +
                        (j - 4 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BBB    = (i + num_ghosts_0_g) +
                        (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_BB     = (i + num_ghosts_0_g) +
                        (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_B      = (i + num_ghosts_0_g) +
                        (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_T      = (i + num_ghosts_0_g) +
                        (j + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TT     = (i + num_ghosts_0_g) +
                        (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTT    = (i + num_ghosts_0_g) +
                        (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTTT   = (i + num_ghosts_0_g) +
                        (j + 3 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTTTT  = (i + num_ghosts_0_g) +
                        (j + 4 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_g_TTTTTT = (i + num_ghosts_0_g) +
                        (j + 5 + num_ghosts_1_g)*ghostcell_dim_0_g;
                    
                    const int idx_h_BBBBBB = (i + num_ghosts_0_h) +
                        (j - 6 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_BBBBB  = (i + num_ghosts_0_h) +
                        (j - 5 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_BBBB   = (i + num_ghosts_0_h) +
                        (j - 4 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_BBB    = (i + num_ghosts_0_h) +
                        (j - 3 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_BB     = (i + num_ghosts_0_h) +
                        (j - 2 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_B      = (i + num_ghosts_0_h) +
                        (j - 1 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_T      = (i + num_ghosts_0_h) +
                        (j + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TT     = (i + num_ghosts_0_h) +
                        (j + 1 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TTT    = (i + num_ghosts_0_h) +
                        (j + 2 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TTTT   = (i + num_ghosts_0_h) +
                        (j + 3 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TTTTT  = (i + num_ghosts_0_h) +
                        (j + 4 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    const int idx_h_TTTTTT = (i + num_ghosts_0_h) +
                        (j + 5 + num_ghosts_1_h)*ghostcell_dim_0_h;
                    
                    F_face_y[idx_face_y] += dt*quarter*(
                        d_coef_a*((f[idx_f_B]      + f[idx_f_T]     )*(g[idx_g_B]      + g[idx_g_T]     )*(h[idx_h_B]      + h[idx_h_T]     )) +
                        d_coef_b*((f[idx_f_BB]     + f[idx_f_T]     )*(g[idx_g_BB]     + g[idx_g_T]     )*(h[idx_h_BB]     + h[idx_h_T]     )  +
                                  (f[idx_f_B]      + f[idx_f_TT]    )*(g[idx_g_B]      + g[idx_g_TT]    )*(h[idx_h_B]      + h[idx_h_TT]    )) +
                        d_coef_c*((f[idx_f_BBB]    + f[idx_f_T]     )*(g[idx_g_BBB]    + g[idx_g_T]     )*(h[idx_h_BBB]    + h[idx_h_T]     )  +
                                  (f[idx_f_BB]     + f[idx_f_TT]    )*(g[idx_g_BB]     + g[idx_g_TT]    )*(h[idx_h_BB]     + h[idx_h_TT]    )  +
                                  (f[idx_f_B]      + f[idx_f_TTT]   )*(g[idx_g_B]      + g[idx_g_TTT]   )*(h[idx_h_B]      + h[idx_h_TTT]   )) +
                        d_coef_d*((f[idx_f_BBBB]   + f[idx_f_T]     )*(g[idx_g_BBBB]   + g[idx_g_T]     )*(h[idx_h_BBBB]   + h[idx_h_T]     )  +
                                  (f[idx_f_BBB]    + f[idx_f_TT]    )*(g[idx_g_BBB]    + g[idx_g_TT]    )*(h[idx_h_BBB]    + h[idx_h_TT]    )  +
                                  (f[idx_f_BB]     + f[idx_f_TTT]   )*(g[idx_g_BB]     + g[idx_g_TTT]   )*(h[idx_h_BB]     + h[idx_h_TTT]   )  +
                                  (f[idx_f_B]      + f[idx_f_TTTT]  )*(g[idx_g_B]      + g[idx_g_TTTT]  )*(h[idx_h_B]      + h[idx_h_TTTT]  )) +
                        d_coef_e*((f[idx_f_BBBBB]  + f[idx_f_T]     )*(g[idx_g_BBBBB]  + g[idx_g_T]     )*(h[idx_h_BBBBB]  + h[idx_h_T]     )  +
                                  (f[idx_f_BBBB]   + f[idx_f_TT]    )*(g[idx_g_BBBB]   + g[idx_g_TT]    )*(h[idx_h_BBBB]   + h[idx_h_TT]    )  +
                                  (f[idx_f_BBB]    + f[idx_f_TTT]   )*(g[idx_g_BBB]    + g[idx_g_TTT]   )*(h[idx_h_BBB]    + h[idx_h_TTT]   )  +
                                  (f[idx_f_BB]     + f[idx_f_TTTT]  )*(g[idx_g_BB]     + g[idx_g_TTTT]  )*(h[idx_h_BB]     + h[idx_h_TTTT]  )  +
                                  (f[idx_f_B]      + f[idx_f_TTTTT] )*(g[idx_g_B]      + g[idx_g_TTTTT] )*(h[idx_h_B]      + h[idx_h_TTTTT] )) +
                        d_coef_f*((f[idx_f_BBBBBB] + f[idx_f_T]     )*(g[idx_g_BBBBBB] + g[idx_g_T]     )*(h[idx_h_BBBBBB] + h[idx_h_T]     )  +
                                  (f[idx_f_BBBBB]  + f[idx_f_TT]    )*(g[idx_g_BBBBB]  + g[idx_g_TT]    )*(h[idx_h_BBBBB]  + h[idx_h_TT]    )  +
                                  (f[idx_f_BBBB]   + f[idx_f_TTT]   )*(g[idx_g_BBBB]   + g[idx_g_TTT]   )*(h[idx_h_BBBB]   + h[idx_h_TTT]   )  +
                                  (f[idx_f_BBB]    + f[idx_f_TTTT]  )*(g[idx_g_BBB]    + g[idx_g_TTTT]  )*(h[idx_h_BBB]    + h[idx_h_TTTT]  )  +
                                  (f[idx_f_BB]     + f[idx_f_TTTTT] )*(g[idx_g_BB]     + g[idx_g_TTTTT] )*(h[idx_h_BB]     + h[idx_h_TTTTT] )  +
                                  (f[idx_f_B]      + f[idx_f_TTTTTT])*(g[idx_g_B]      + g[idx_g_TTTTTT])*(h[idx_h_B]      + h[idx_h_TTTTTT])));
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
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_B = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_B = (i + num_ghosts_0_g) +
                            (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_T = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_B = (i + num_ghosts_0_h) +
                            (j - 1 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_T = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_y[idx_face_y] += dt*quarter*(
                            d_coef_a*(f[idx_f_B] + f[idx_f_T])*(g[idx_g_B] + g[idx_g_T])*(h[idx_h_B] + h[idx_h_T]));
                    }
                }
            }
        }
        else if (d_stencil_width == 5)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BB = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B  = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BB = (i + num_ghosts_0_g) +
                            (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B  = (i + num_ghosts_0_g) +
                            (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_T  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TT = (i + num_ghosts_0_g) +
                            (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_BB = (i + num_ghosts_0_h) +
                            (j - 2 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_B  = (i + num_ghosts_0_h) +
                            (j - 1 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_T  = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TT = (i + num_ghosts_0_h) +
                            (j + 1 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_y[idx_face_y] += dt*quarter*(
                            d_coef_a*((f[idx_f_B]  + f[idx_f_T] )*(g[idx_g_B]  + g[idx_g_T] )*(h[idx_h_B]  + h[idx_h_T] )) +
                            d_coef_b*((f[idx_f_BB] + f[idx_f_T] )*(g[idx_g_BB] + g[idx_g_T] )*(h[idx_h_BB] + h[idx_h_T] )  +
                                      (f[idx_f_B]  + f[idx_f_TT])*(g[idx_g_B]  + g[idx_g_TT])*(h[idx_h_B]  + h[idx_h_TT])));
                    }
                }
            }
        }
        else if (d_stencil_width == 7)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BBB = (i + num_ghosts_0_f) +
                            (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB  = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B   = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT  = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTT = (i + num_ghosts_0_f) +
                            (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBB = (i + num_ghosts_0_g) +
                            (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB  = (i + num_ghosts_0_g) +
                            (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B   = (i + num_ghosts_0_g) +
                            (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_T   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TT  = (i + num_ghosts_0_g) +
                            (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTT = (i + num_ghosts_0_g) +
                            (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_BBB = (i + num_ghosts_0_h) +
                            (j - 3 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BB  = (i + num_ghosts_0_h) +
                            (j - 2 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_B   = (i + num_ghosts_0_h) +
                            (j - 1 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_T   = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TT  = (i + num_ghosts_0_h) +
                            (j + 1 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TTT = (i + num_ghosts_0_h) +
                            (j + 2 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_y[idx_face_y] += dt*quarter*(
                            d_coef_a*((f[idx_f_B]   + f[idx_f_T]  )*(g[idx_g_B]   + g[idx_g_T]  )*(h[idx_h_B]   + h[idx_h_T]  )) +
                            d_coef_b*((f[idx_f_BB]  + f[idx_f_T]  )*(g[idx_g_BB]  + g[idx_g_T]  )*(h[idx_h_BB]  + h[idx_h_T]  )  +
                                      (f[idx_f_B]   + f[idx_f_TT] )*(g[idx_g_B]   + g[idx_g_TT] )*(h[idx_h_B]   + h[idx_h_TT] )) +
                            d_coef_c*((f[idx_f_BBB] + f[idx_f_T]  )*(g[idx_g_BBB] + g[idx_g_T]  )*(h[idx_h_BBB] + h[idx_h_T]  )  +
                                      (f[idx_f_BB]  + f[idx_f_TT] )*(g[idx_g_BB]  + g[idx_g_TT] )*(h[idx_h_BB]  + h[idx_h_TT] )  +
                                      (f[idx_f_B]   + f[idx_f_TTT])*(g[idx_g_B]   + g[idx_g_TTT])*(h[idx_h_B]   + h[idx_h_TTT])));
                    }
                }
            }
        }
        else if (d_stencil_width == 9)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BBBB = (i + num_ghosts_0_f) +
                            (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB  = (i + num_ghosts_0_f) +
                            (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB   = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B    = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT   = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTT  = (i + num_ghosts_0_f) +
                            (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTT = (i + num_ghosts_0_f) +
                            (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBBB = (i + num_ghosts_0_g) +
                            (j - 4 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBB  = (i + num_ghosts_0_g) +
                            (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB   = (i + num_ghosts_0_g) +
                            (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B    = (i + num_ghosts_0_g) +
                            (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_T    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TT   = (i + num_ghosts_0_g) +
                            (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTT  = (i + num_ghosts_0_g) +
                            (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTTT = (i + num_ghosts_0_g) +
                            (j + 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_BBBB = (i + num_ghosts_0_h) +
                            (j - 4 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BBB  = (i + num_ghosts_0_h) +
                            (j - 3 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BB   = (i + num_ghosts_0_h) +
                            (j - 2 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_B    = (i + num_ghosts_0_h) +
                            (j - 1 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_T    = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TT   = (i + num_ghosts_0_h) +
                            (j + 1 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TTT  = (i + num_ghosts_0_h) +
                            (j + 2 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TTTT = (i + num_ghosts_0_h) +
                            (j + 3 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_y[idx_face_y] += dt*quarter*(
                            d_coef_a*((f[idx_f_B]    + f[idx_f_T]   )*(g[idx_g_B]    + g[idx_g_T]   )*(h[idx_h_B]    + h[idx_h_T]   )) +
                            d_coef_b*((f[idx_f_BB]   + f[idx_f_T]   )*(g[idx_g_BB]   + g[idx_g_T]   )*(h[idx_h_BB]   + h[idx_h_T]   )  +
                                      (f[idx_f_B]    + f[idx_f_TT]  )*(g[idx_g_B]    + g[idx_g_TT]  )*(h[idx_h_B]    + h[idx_h_TT]  )) +
                            d_coef_c*((f[idx_f_BBB]  + f[idx_f_T]   )*(g[idx_g_BBB]  + g[idx_g_T]   )*(h[idx_h_BBB]  + h[idx_h_T]   )  +
                                      (f[idx_f_BB]   + f[idx_f_TT]  )*(g[idx_g_BB]   + g[idx_g_TT]  )*(h[idx_h_BB]   + h[idx_h_TT]  )  +
                                      (f[idx_f_B]    + f[idx_f_TTT] )*(g[idx_g_B]    + g[idx_g_TTT] )*(h[idx_h_B]    + h[idx_h_TTT] )) +
                            d_coef_d*((f[idx_f_BBBB] + f[idx_f_T]   )*(g[idx_g_BBBB] + g[idx_g_T]   )*(h[idx_h_BBBB] + h[idx_h_T]   )  +
                                      (f[idx_f_BBB]  + f[idx_f_TT]  )*(g[idx_g_BBB]  + g[idx_g_TT]  )*(h[idx_h_BBB]  + h[idx_h_TT]  )  +
                                      (f[idx_f_BB]   + f[idx_f_TTT] )*(g[idx_g_BB]   + g[idx_g_TTT] )*(h[idx_h_BB]   + h[idx_h_TTT] )  +
                                      (f[idx_f_B]    + f[idx_f_TTTT])*(g[idx_g_B]    + g[idx_g_TTTT])*(h[idx_h_B]    + h[idx_h_TTTT])));
                    }
                }
            }
        }
        else if (d_stencil_width == 11)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BBBBB = (i + num_ghosts_0_f) +
                            (j - 5 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBB  = (i + num_ghosts_0_f) +
                            (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB   = (i + num_ghosts_0_f) +
                            (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB    = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B     = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT    = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTT   = (i + num_ghosts_0_f) +
                            (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTT  = (i + num_ghosts_0_f) +
                            (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTTT = (i + num_ghosts_0_f) +
                            (j + 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBBBB = (i + num_ghosts_0_g) +
                            (j - 5 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBBB  = (i + num_ghosts_0_g) +
                            (j - 4 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBB   = (i + num_ghosts_0_g) +
                            (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB    = (i + num_ghosts_0_g) +
                            (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B     = (i + num_ghosts_0_g) +
                            (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_T     = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TT    = (i + num_ghosts_0_g) +
                            (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTT   = (i + num_ghosts_0_g) +
                            (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTTT  = (i + num_ghosts_0_g) +
                            (j + 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTTTT = (i + num_ghosts_0_g) +
                            (j + 4 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_BBBBB = (i + num_ghosts_0_h) +
                            (j - 5 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BBBB  = (i + num_ghosts_0_h) +
                            (j - 4 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BBB   = (i + num_ghosts_0_h) +
                            (j - 3 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BB    = (i + num_ghosts_0_h) +
                            (j - 2 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_B     = (i + num_ghosts_0_h) +
                            (j - 1 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_T     = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TT    = (i + num_ghosts_0_h) +
                            (j + 1 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TTT   = (i + num_ghosts_0_h) +
                            (j + 2 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TTTT  = (i + num_ghosts_0_h) +
                            (j + 3 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TTTTT = (i + num_ghosts_0_h) +
                            (j + 4 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_y[idx_face_y] += dt*quarter*(
                            d_coef_a*((f[idx_f_B]     + f[idx_f_T]    )*(g[idx_g_B]     + g[idx_g_T]    )*(h[idx_h_B]     + h[idx_h_T]    )) +
                            d_coef_b*((f[idx_f_BB]    + f[idx_f_T]    )*(g[idx_g_BB]    + g[idx_g_T]    )*(h[idx_h_BB]    + h[idx_h_T]    )  +
                                      (f[idx_f_B]     + f[idx_f_TT]   )*(g[idx_g_B]     + g[idx_g_TT]   )*(h[idx_h_B]     + h[idx_h_TT]   )) +
                            d_coef_c*((f[idx_f_BBB]   + f[idx_f_T]    )*(g[idx_g_BBB]   + g[idx_g_T]    )*(h[idx_h_BBB]   + h[idx_h_T]    )  +
                                      (f[idx_f_BB]    + f[idx_f_TT]   )*(g[idx_g_BB]    + g[idx_g_TT]   )*(h[idx_h_BB]    + h[idx_h_TT]   )  +
                                      (f[idx_f_B]     + f[idx_f_TTT]  )*(g[idx_g_B]     + g[idx_g_TTT]  )*(h[idx_h_B]     + h[idx_h_TTT]  )) +
                            d_coef_d*((f[idx_f_BBBB]  + f[idx_f_T]    )*(g[idx_g_BBBB]  + g[idx_g_T]    )*(h[idx_h_BBBB]  + h[idx_h_T]    )  +
                                      (f[idx_f_BBB]   + f[idx_f_TT]   )*(g[idx_g_BBB]   + g[idx_g_TT]   )*(h[idx_h_BBB]   + h[idx_h_TT]   )  +
                                      (f[idx_f_BB]    + f[idx_f_TTT]  )*(g[idx_g_BB]    + g[idx_g_TTT]  )*(h[idx_h_BB]    + h[idx_h_TTT]  )  +
                                      (f[idx_f_B]     + f[idx_f_TTTT] )*(g[idx_g_B]     + g[idx_g_TTTT] )*(h[idx_h_B]     + h[idx_h_TTTT] )) +
                            d_coef_e*((f[idx_f_BBBBB] + f[idx_f_T]    )*(g[idx_g_BBBBB] + g[idx_g_T]    )*(h[idx_h_BBBBB] + h[idx_h_T]    )  +
                                      (f[idx_f_BBBB]  + f[idx_f_TT]   )*(g[idx_g_BBBB]  + g[idx_g_TT]   )*(h[idx_h_BBBB]  + h[idx_h_TT]   )  +
                                      (f[idx_f_BBB]   + f[idx_f_TTT]  )*(g[idx_g_BBB]   + g[idx_g_TTT]  )*(h[idx_h_BBB]   + h[idx_h_TTT]  )  +
                                      (f[idx_f_BB]    + f[idx_f_TTTT] )*(g[idx_g_BB]    + g[idx_g_TTTT] )*(h[idx_h_BB]    + h[idx_h_TTTT] )  +
                                      (f[idx_f_B]     + f[idx_f_TTTTT])*(g[idx_g_B]     + g[idx_g_TTTTT])*(h[idx_h_B]     + h[idx_h_TTTTT])));
                    }
                }
            }
        }
        else if (d_stencil_width == 13)
        {
            for (int k = 0; k < interior_dim_2; k++)
            {
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_y = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*(interior_dim_1 + 1);
                        
                        const int idx_f_BBBBBB = (i + num_ghosts_0_f) +
                            (j - 6 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBBB  = (i + num_ghosts_0_f) +
                            (j - 5 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBB   = (i + num_ghosts_0_f) +
                            (j - 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB    = (i + num_ghosts_0_f) +
                            (j - 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB     = (i + num_ghosts_0_f) +
                            (j - 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B      = (i + num_ghosts_0_f) +
                            (j - 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_T      = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TT     = (i + num_ghosts_0_f) +
                            (j + 1 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTT    = (i + num_ghosts_0_f) +
                            (j + 2 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTT   = (i + num_ghosts_0_f) +
                            (j + 3 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTTT  = (i + num_ghosts_0_f) +
                            (j + 4 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_TTTTTT = (i + num_ghosts_0_f) +
                            (j + 5 + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBBBBB = (i + num_ghosts_0_g) +
                            (j - 6 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBBBB  = (i + num_ghosts_0_g) +
                            (j - 5 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBBB   = (i + num_ghosts_0_g) +
                            (j - 4 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBB    = (i + num_ghosts_0_g) +
                            (j - 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB     = (i + num_ghosts_0_g) +
                            (j - 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B      = (i + num_ghosts_0_g) +
                            (j - 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_T      = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TT     = (i + num_ghosts_0_g) +
                            (j + 1 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTT    = (i + num_ghosts_0_g) +
                            (j + 2 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTTT   = (i + num_ghosts_0_g) +
                            (j + 3 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTTTT  = (i + num_ghosts_0_g) +
                            (j + 4 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_TTTTTT = (i + num_ghosts_0_g) +
                            (j + 5 + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_BBBBBB = (i + num_ghosts_0_h) +
                            (j - 6 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BBBBB  = (i + num_ghosts_0_h) +
                            (j - 5 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BBBB   = (i + num_ghosts_0_h) +
                            (j - 4 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BBB    = (i + num_ghosts_0_h) +
                            (j - 3 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BB     = (i + num_ghosts_0_h) +
                            (j - 2 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_B      = (i + num_ghosts_0_h) +
                            (j - 1 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_T      = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TT     = (i + num_ghosts_0_h) +
                            (j + 1 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TTT    = (i + num_ghosts_0_h) +
                            (j + 2 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TTTT   = (i + num_ghosts_0_h) +
                            (j + 3 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TTTTT  = (i + num_ghosts_0_h) +
                            (j + 4 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_TTTTTT = (i + num_ghosts_0_h) +
                            (j + 5 + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_y[idx_face_y] += dt*quarter*(
                            d_coef_a*((f[idx_f_B]      + f[idx_f_T]     )*(g[idx_g_B]      + g[idx_g_T]     )*(h[idx_h_B]      + h[idx_h_T]     )) +
                            d_coef_b*((f[idx_f_BB]     + f[idx_f_T]     )*(g[idx_g_BB]     + g[idx_g_T]     )*(h[idx_h_BB]     + h[idx_h_T]     )  +
                                      (f[idx_f_B]      + f[idx_f_TT]    )*(g[idx_g_B]      + g[idx_g_TT]    )*(h[idx_h_B]      + h[idx_h_TT]    )) +
                            d_coef_c*((f[idx_f_BBB]    + f[idx_f_T]     )*(g[idx_g_BBB]    + g[idx_g_T]     )*(h[idx_h_BBB]    + h[idx_h_T]     )  +
                                      (f[idx_f_BB]     + f[idx_f_TT]    )*(g[idx_g_BB]     + g[idx_g_TT]    )*(h[idx_h_BB]     + h[idx_h_TT]    )  +
                                      (f[idx_f_B]      + f[idx_f_TTT]   )*(g[idx_g_B]      + g[idx_g_TTT]   )*(h[idx_h_B]      + h[idx_h_TTT]   )) +
                            d_coef_d*((f[idx_f_BBBB]   + f[idx_f_T]     )*(g[idx_g_BBBB]   + g[idx_g_T]     )*(h[idx_h_BBBB]   + h[idx_h_T]     )  +
                                      (f[idx_f_BBB]    + f[idx_f_TT]    )*(g[idx_g_BBB]    + g[idx_g_TT]    )*(h[idx_h_BBB]    + h[idx_h_TT]    )  +
                                      (f[idx_f_BB]     + f[idx_f_TTT]   )*(g[idx_g_BB]     + g[idx_g_TTT]   )*(h[idx_h_BB]     + h[idx_h_TTT]   )  +
                                      (f[idx_f_B]      + f[idx_f_TTTT]  )*(g[idx_g_B]      + g[idx_g_TTTT]  )*(h[idx_h_B]      + h[idx_h_TTTT]  )) +
                            d_coef_e*((f[idx_f_BBBBB]  + f[idx_f_T]     )*(g[idx_g_BBBBB]  + g[idx_g_T]     )*(h[idx_h_BBBBB]  + h[idx_h_T]     )  +
                                      (f[idx_f_BBBB]   + f[idx_f_TT]    )*(g[idx_g_BBBB]   + g[idx_g_TT]    )*(h[idx_h_BBBB]   + h[idx_h_TT]    )  +
                                      (f[idx_f_BBB]    + f[idx_f_TTT]   )*(g[idx_g_BBB]    + g[idx_g_TTT]   )*(h[idx_h_BBB]    + h[idx_h_TTT]   )  +
                                      (f[idx_f_BB]     + f[idx_f_TTTT]  )*(g[idx_g_BB]     + g[idx_g_TTTT]  )*(h[idx_h_BB]     + h[idx_h_TTTT]  )  +
                                      (f[idx_f_B]      + f[idx_f_TTTTT] )*(g[idx_g_B]      + g[idx_g_TTTTT] )*(h[idx_h_B]      + h[idx_h_TTTTT] )) +
                            d_coef_f*((f[idx_f_BBBBBB] + f[idx_f_T]     )*(g[idx_g_BBBBBB] + g[idx_g_T]     )*(h[idx_h_BBBBBB] + h[idx_h_T]     )  +
                                      (f[idx_f_BBBBB]  + f[idx_f_TT]    )*(g[idx_g_BBBBB]  + g[idx_g_TT]    )*(h[idx_h_BBBBB]  + h[idx_h_TT]    )  +
                                      (f[idx_f_BBBB]   + f[idx_f_TTT]   )*(g[idx_g_BBBB]   + g[idx_g_TTT]   )*(h[idx_h_BBBB]   + h[idx_h_TTT]   )  +
                                      (f[idx_f_BBB]    + f[idx_f_TTTT]  )*(g[idx_g_BBB]    + g[idx_g_TTTT]  )*(h[idx_h_BBB]    + h[idx_h_TTTT]  )  +
                                      (f[idx_f_BB]     + f[idx_f_TTTTT] )*(g[idx_g_BB]     + g[idx_g_TTTTT] )*(h[idx_h_BB]     + h[idx_h_TTTTT] )  +
                                      (f[idx_f_B]      + f[idx_f_TTTTTT])*(g[idx_g_B]      + g[idx_g_TTTTTT])*(h[idx_h_B]      + h[idx_h_TTTTTT])));
                    }
                }
            }
        }
    }
}


/*
 * Add linear term to convective flux in z-direction.
 */
void
ConvectiveFluxReconstructorKEP::addLinearTermToConvectiveFluxZ(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_f,
    const int component_idx_flux,
    const int component_idx_f,
    const double dt) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR("ConvectiveFluxReconstructorKEP::addLinearTermToConvectiveFluxZ:"
            " Cannot compute linear term in z-direction for 1D problem!");
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR("ConvectiveFluxReconstructorKEP::addLinearTermToConvectiveFluxZ:"
            " Cannot compute linear term in z-direction for 2D problem!");
    }
#endif
    
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
    
    double* F_face_z = data_convective_flux->getPointer(2, component_idx_flux);
    double* f        = data_f->getPointer(component_idx_f);
    
    if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        /*
         * Get the numbers of ghost cells and the dimensions of the ghost cell box.
         */
        
        const int num_ghosts_0_f = num_ghosts_f[0];
        const int num_ghosts_1_f = num_ghosts_f[1];
        const int num_ghosts_2_f = num_ghosts_f[2];
        const int ghostcell_dim_0_f = ghostcell_dims_f[0];
        const int ghostcell_dim_1_f = ghostcell_dims_f[1];
        
        
        if (d_stencil_width == 3)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_B = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_z[idx_face_z] += dt*(
                            d_coef_a*(f[idx_f_B] + f[idx_f_F]));
                    }
                }
            }
        }
        else if (d_stencil_width == 5)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_z[idx_face_z] += dt*(
                            d_coef_a*((f[idx_f_B]  + f[idx_f_F] )) +
                            d_coef_b*((f[idx_f_BB] + f[idx_f_F] )  +
                                      (f[idx_f_B]  + f[idx_f_FF])));
                    }
                }
            }
        }
        else if (d_stencil_width == 7)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BBB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_z[idx_face_z] += dt*(
                            d_coef_a*((f[idx_f_B]   + f[idx_f_F]  )) +
                            d_coef_b*((f[idx_f_BB]  + f[idx_f_F]  )  +
                                      (f[idx_f_B]   + f[idx_f_FF] )) +
                            d_coef_c*((f[idx_f_BBB] + f[idx_f_F]  )  +
                                      (f[idx_f_BB]  + f[idx_f_FF] )  +
                                      (f[idx_f_B]   + f[idx_f_FFF])));
                    }
                }
            }
        }
        else if (d_stencil_width == 9)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BBBB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFF  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_z[idx_face_z] += dt*(
                            d_coef_a*((f[idx_f_B]    + f[idx_f_F]   )) +
                            d_coef_b*((f[idx_f_BB]   + f[idx_f_F]   )  +
                                      (f[idx_f_B]    + f[idx_f_FF]  )) +
                            d_coef_c*((f[idx_f_BBB]  + f[idx_f_F]   )  +
                                      (f[idx_f_BB]   + f[idx_f_FF]  )  +
                                      (f[idx_f_B]    + f[idx_f_FFF] )) +
                            d_coef_d*((f[idx_f_BBBB] + f[idx_f_F]   )  +
                                      (f[idx_f_BBB]  + f[idx_f_FF]  )  +
                                      (f[idx_f_BB]   + f[idx_f_FFF] )  +
                                      (f[idx_f_B]    + f[idx_f_FFFF])));
                    }
                }
            }
        }
        else if (d_stencil_width == 11)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BBBBB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 5 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBB  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFF   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFF  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFFF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_z[idx_face_z] += dt*(
                            d_coef_a*((f[idx_f_B]     + f[idx_f_F]    )) +
                            d_coef_b*((f[idx_f_BB]    + f[idx_f_F]    )  +
                                      (f[idx_f_B]     + f[idx_f_FF]   )) +
                            d_coef_c*((f[idx_f_BBB]   + f[idx_f_F]    )  +
                                      (f[idx_f_BB]    + f[idx_f_FF]   )  +
                                      (f[idx_f_B]     + f[idx_f_FFF]  )) +
                            d_coef_d*((f[idx_f_BBBB]  + f[idx_f_F]    )  +
                                      (f[idx_f_BBB]   + f[idx_f_FF]   )  +
                                      (f[idx_f_BB]    + f[idx_f_FFF]  )  +
                                      (f[idx_f_B]     + f[idx_f_FFFF] )) +
                            d_coef_e*((f[idx_f_BBBBB] + f[idx_f_F]    )  +
                                      (f[idx_f_BBBB]  + f[idx_f_FF]   )  +
                                      (f[idx_f_BBB]   + f[idx_f_FFF]  )  +
                                      (f[idx_f_BB]    + f[idx_f_FFFF] )  +
                                      (f[idx_f_B]     + f[idx_f_FFFFF])));
                    }
                }
            }
        }
        else if (d_stencil_width == 13)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BBBBBB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 6 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBBB  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 5 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBB   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B      = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F      = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFF    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFF   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFFF  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFFFF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 5 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        F_face_z[idx_face_z] += dt*(
                            d_coef_a*((f[idx_f_B]      + f[idx_f_F]     )) +
                            d_coef_b*((f[idx_f_BB]     + f[idx_f_F]     )  +
                                      (f[idx_f_B]      + f[idx_f_FF]    )) +
                            d_coef_c*((f[idx_f_BBB]    + f[idx_f_F]     )  +
                                      (f[idx_f_BB]     + f[idx_f_FF]    )  +
                                      (f[idx_f_B]      + f[idx_f_FFF]   )) +
                            d_coef_d*((f[idx_f_BBBB]   + f[idx_f_F]     )  +
                                      (f[idx_f_BBB]    + f[idx_f_FF]    )  +
                                      (f[idx_f_BB]     + f[idx_f_FFF]   )  +
                                      (f[idx_f_B]      + f[idx_f_FFFF]  )) +
                            d_coef_e*((f[idx_f_BBBBB]  + f[idx_f_F]     )  +
                                      (f[idx_f_BBBB]   + f[idx_f_FF]    )  +
                                      (f[idx_f_BBB]    + f[idx_f_FFF]   )  +
                                      (f[idx_f_BB]     + f[idx_f_FFFF]  )  +
                                      (f[idx_f_B]      + f[idx_f_FFFFF] )) +
                            d_coef_f*((f[idx_f_BBBBBB] + f[idx_f_F]     )  +
                                      (f[idx_f_BBBBB]  + f[idx_f_FF]    )  +
                                      (f[idx_f_BBBB]   + f[idx_f_FFF]   )  +
                                      (f[idx_f_BBB]    + f[idx_f_FFFF]  )  +
                                      (f[idx_f_BB]     + f[idx_f_FFFFF] )  +
                                      (f[idx_f_B]      + f[idx_f_FFFFFF])));
                    }
                }
            }
        }
    }
}


/*
 * Add quadratic term to convective flux in z-direction.
 */
void
ConvectiveFluxReconstructorKEP::addQuadraticTermToConvectiveFluxZ(
    HAMERS_SHARED_PTR<pdat::SideData<double> >& data_convective_flux,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_f,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& data_g,
    const int component_idx_flux,
    const int component_idx_f,
    const int component_idx_g,
    const double dt) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR("ConvectiveFluxReconstructorKEP::addQuadraticTermToConvectiveFluxZ:"
            " Cannot compute quadratic term in z-direction for 1D problem!");
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR("ConvectiveFluxReconstructorKEP::addQuadraticTermToConvectiveFluxZ:"
            " Cannot compute quadratic term in z-direction for 2D problem!");
    }
#endif
    
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
    
    /*
     * Get the pointers to the data.
     */
    
    double* F_face_z = data_convective_flux->getPointer(2, component_idx_flux);
    double* f        = data_f->getPointer(component_idx_f);
    double* g        = data_g->getPointer(component_idx_g);
    
    const double half = double(1)/double(2);
    
    if (d_dim == tbox::Dimension(3))
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
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_B = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_B = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_F = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_z[idx_face_z] += dt*half*(
                            d_coef_a*(f[idx_f_B] + f[idx_f_F])*(g[idx_g_B] + g[idx_g_F]));
                    }
                }
            }
        }
        else if (d_stencil_width == 5)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BB = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_F  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FF = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_z[idx_face_z] += dt*half*(
                            d_coef_a*((f[idx_f_B]  + f[idx_f_F] )*(g[idx_g_B]  + g[idx_g_F] )) +
                            d_coef_b*((f[idx_f_BB] + f[idx_f_F] )*(g[idx_g_BB] + g[idx_g_F] )  +
                                      (f[idx_f_B]  + f[idx_f_FF])*(g[idx_g_B]  + g[idx_g_FF])));
                    }
                }
            }
        }
        else if (d_stencil_width == 7)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BBB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBB = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_F   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FF  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFF = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_z[idx_face_z] += dt*half*(
                            d_coef_a*((f[idx_f_B]   + f[idx_f_F]  )*(g[idx_g_B]   + g[idx_g_F]  )) +
                            d_coef_b*((f[idx_f_BB]  + f[idx_f_F]  )*(g[idx_g_BB]  + g[idx_g_F]  )  +
                                      (f[idx_f_B]   + f[idx_f_FF] )*(g[idx_g_B]   + g[idx_g_FF] )) +
                            d_coef_c*((f[idx_f_BBB] + f[idx_f_F]  )*(g[idx_g_BBB] + g[idx_g_F]  )  +
                                      (f[idx_f_BB]  + f[idx_f_FF] )*(g[idx_g_BB]  + g[idx_g_FF] )  +
                                      (f[idx_f_B]   + f[idx_f_FFF])*(g[idx_g_B]   + g[idx_g_FFF])));
                    }
                }
            }
        }
        else if (d_stencil_width == 9)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BBBB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFF  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBBB = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 4 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBB  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_F    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FF   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFF  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFFF = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_z[idx_face_z] += dt*half*(
                            d_coef_a*((f[idx_f_B]    + f[idx_f_F]   )*(g[idx_g_B]    + g[idx_g_F]   )) +
                            d_coef_b*((f[idx_f_BB]   + f[idx_f_F]   )*(g[idx_g_BB]   + g[idx_g_F]   )  +
                                      (f[idx_f_B]    + f[idx_f_FF]  )*(g[idx_g_B]    + g[idx_g_FF]  )) +
                            d_coef_c*((f[idx_f_BBB]  + f[idx_f_F]   )*(g[idx_g_BBB]  + g[idx_g_F]   )  +
                                      (f[idx_f_BB]   + f[idx_f_FF]  )*(g[idx_g_BB]   + g[idx_g_FF]  )  +
                                      (f[idx_f_B]    + f[idx_f_FFF] )*(g[idx_g_B]    + g[idx_g_FFF] )) +
                            d_coef_d*((f[idx_f_BBBB] + f[idx_f_F]   )*(g[idx_g_BBBB] + g[idx_g_F]   )  +
                                      (f[idx_f_BBB]  + f[idx_f_FF]  )*(g[idx_g_BBB]  + g[idx_g_FF]  )  +
                                      (f[idx_f_BB]   + f[idx_f_FFF] )*(g[idx_g_BB]   + g[idx_g_FFF] )  +
                                      (f[idx_f_B]    + f[idx_f_FFFF])*(g[idx_g_B]    + g[idx_g_FFFF])));
                    }
                }
            }
        }
        else if (d_stencil_width == 11)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BBBBB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 5 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBB  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFF   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFF  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFFF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBBBB = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 5 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBBB  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 4 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBB   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B     = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_F     = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FF    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFF   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFFF  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFFFF = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 4 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_z[idx_face_z] += dt*half*(
                            d_coef_a*((f[idx_f_B]     + f[idx_f_F]    )*(g[idx_g_B]     + g[idx_g_F]    )) +
                            d_coef_b*((f[idx_f_BB]    + f[idx_f_F]    )*(g[idx_g_BB]    + g[idx_g_F]    )  +
                                      (f[idx_f_B]     + f[idx_f_FF]   )*(g[idx_g_B]     + g[idx_g_FF]   )) +
                            d_coef_c*((f[idx_f_BBB]   + f[idx_f_F]    )*(g[idx_g_BBB]   + g[idx_g_F]    )  +
                                      (f[idx_f_BB]    + f[idx_f_FF]   )*(g[idx_g_BB]    + g[idx_g_FF]   )  +
                                      (f[idx_f_B]     + f[idx_f_FFF]  )*(g[idx_g_B]     + g[idx_g_FFF]  )) +
                            d_coef_d*((f[idx_f_BBBB]  + f[idx_f_F]    )*(g[idx_g_BBBB]  + g[idx_g_F]    )  +
                                      (f[idx_f_BBB]   + f[idx_f_FF]   )*(g[idx_g_BBB]   + g[idx_g_FF]   )  +
                                      (f[idx_f_BB]    + f[idx_f_FFF]  )*(g[idx_g_BB]    + g[idx_g_FFF]  )  +
                                      (f[idx_f_B]     + f[idx_f_FFFF] )*(g[idx_g_B]     + g[idx_g_FFFF] )) +
                            d_coef_e*((f[idx_f_BBBBB] + f[idx_f_F]    )*(g[idx_g_BBBBB] + g[idx_g_F]    )  +
                                      (f[idx_f_BBBB]  + f[idx_f_FF]   )*(g[idx_g_BBBB]  + g[idx_g_FF]   )  +
                                      (f[idx_f_BBB]   + f[idx_f_FFF]  )*(g[idx_g_BBB]   + g[idx_g_FFF]  )  +
                                      (f[idx_f_BB]    + f[idx_f_FFFF] )*(g[idx_g_BB]    + g[idx_g_FFFF] )  +
                                      (f[idx_f_B]     + f[idx_f_FFFFF])*(g[idx_g_B]     + g[idx_g_FFFFF])));
                    }
                }
            }
        }
        else if (d_stencil_width == 13)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BBBBBB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 6 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBBB  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 5 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBB   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B      = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F      = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFF    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFF   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFFF  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFFFF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 5 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBBBBB = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 6 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBBBB  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 5 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBBB   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 4 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBB    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB     = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B      = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_F      = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FF     = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFF    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFFF   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFFFF  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 4 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFFFFF = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 5 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        F_face_z[idx_face_z] += dt*half*(
                            d_coef_a*((f[idx_f_B]      + f[idx_f_F]     )*(g[idx_g_B]      + g[idx_g_F]     )) +
                            d_coef_b*((f[idx_f_BB]     + f[idx_f_F]     )*(g[idx_g_BB]     + g[idx_g_F]     )  +
                                      (f[idx_f_B]      + f[idx_f_FF]    )*(g[idx_g_B]      + g[idx_g_FF]    )) +
                            d_coef_c*((f[idx_f_BBB]    + f[idx_f_F]     )*(g[idx_g_BBB]    + g[idx_g_F]     )  +
                                      (f[idx_f_BB]     + f[idx_f_FF]    )*(g[idx_g_BB]     + g[idx_g_FF]    )  +
                                      (f[idx_f_B]      + f[idx_f_FFF]   )*(g[idx_g_B]      + g[idx_g_FFF]   )) +
                            d_coef_d*((f[idx_f_BBBB]   + f[idx_f_F]     )*(g[idx_g_BBBB]   + g[idx_g_F]     )  +
                                      (f[idx_f_BBB]    + f[idx_f_FF]    )*(g[idx_g_BBB]    + g[idx_g_FF]    )  +
                                      (f[idx_f_BB]     + f[idx_f_FFF]   )*(g[idx_g_BB]     + g[idx_g_FFF]   )  +
                                      (f[idx_f_B]      + f[idx_f_FFFF]  )*(g[idx_g_B]      + g[idx_g_FFFF]  )) +
                            d_coef_e*((f[idx_f_BBBBB]  + f[idx_f_F]     )*(g[idx_g_BBBBB]  + g[idx_g_F]     )  +
                                      (f[idx_f_BBBB]   + f[idx_f_FF]    )*(g[idx_g_BBBB]   + g[idx_g_FF]    )  +
                                      (f[idx_f_BBB]    + f[idx_f_FFF]   )*(g[idx_g_BBB]    + g[idx_g_FFF]   )  +
                                      (f[idx_f_BB]     + f[idx_f_FFFF]  )*(g[idx_g_BB]     + g[idx_g_FFFF]  )  +
                                      (f[idx_f_B]      + f[idx_f_FFFFF] )*(g[idx_g_B]      + g[idx_g_FFFFF] )) +
                            d_coef_f*((f[idx_f_BBBBBB] + f[idx_f_F]     )*(g[idx_g_BBBBBB] + g[idx_g_F]     )  +
                                      (f[idx_f_BBBBB]  + f[idx_f_FF]    )*(g[idx_g_BBBBB]  + g[idx_g_FF]    )  +
                                      (f[idx_f_BBBB]   + f[idx_f_FFF]   )*(g[idx_g_BBBB]   + g[idx_g_FFF]   )  +
                                      (f[idx_f_BBB]    + f[idx_f_FFFF]  )*(g[idx_g_BBB]    + g[idx_g_FFFF]  )  +
                                      (f[idx_f_BB]     + f[idx_f_FFFFF] )*(g[idx_g_BB]     + g[idx_g_FFFFF] )  +
                                      (f[idx_f_B]      + f[idx_f_FFFFFF])*(g[idx_g_B]      + g[idx_g_FFFFFF])));
                    }
                }
            }
        }
    }
}


/*
 * Add cubic term to convective flux in z-direction.
 */
void
ConvectiveFluxReconstructorKEP::addCubicTermToConvectiveFluxZ(
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
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR("ConvectiveFluxReconstructorKEP::addCubicTermToConvectiveFluxZ:"
            " Cannot compute cubic term in z-direction for 1D problem!");
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR("ConvectiveFluxReconstructorKEP::addCubicTermToConvectiveFluxZ:"
            " Cannot compute cubic term in z-direction for 2D problem!");
    }
#endif
    
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
    
    double* F_face_z = data_convective_flux->getPointer(2, component_idx_flux);
    double* f        = data_f->getPointer(component_idx_f);
    double* g        = data_g->getPointer(component_idx_g);
    double* h        = data_h->getPointer(component_idx_h);
    
    const double quarter = double(1)/double(4);
    
    if (d_dim == tbox::Dimension(3))
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
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_B = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_B = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_F = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_B = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 1 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_F = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_z[idx_face_z] += dt*quarter*(
                            d_coef_a*(f[idx_f_B] + f[idx_f_F])*(g[idx_g_B] + g[idx_g_F])*(h[idx_h_B] + h[idx_h_F]));
                    }
                }
            }
        }
        else if (d_stencil_width == 5)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BB = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_F  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FF = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_BB = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 2 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_B  = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 1 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_F  = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FF = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 1 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_z[idx_face_z] += dt*quarter*(
                            d_coef_a*((f[idx_f_B]  + f[idx_f_F] )*(g[idx_g_B]  + g[idx_g_F] )*(h[idx_h_B]  + h[idx_h_F] )) +
                            d_coef_b*((f[idx_f_BB] + f[idx_f_F] )*(g[idx_g_BB] + g[idx_g_F] )*(h[idx_h_BB] + h[idx_h_F] )  +
                                      (f[idx_f_B]  + f[idx_f_FF])*(g[idx_g_B]  + g[idx_g_FF])*(h[idx_h_B]  + h[idx_h_FF])));
                    }
                }
            }
        }
        else if (d_stencil_width == 7)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BBB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBB = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_F   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FF  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFF = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_BBB = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 3 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BB  = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 2 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_B   = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 1 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_F   = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FF  = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 1 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FFF = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 2 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_z[idx_face_z] += dt*quarter*(
                            d_coef_a*((f[idx_f_B]   + f[idx_f_F]  )*(g[idx_g_B]   + g[idx_g_F]  )*(h[idx_h_B]   + h[idx_h_F]  )) +
                            d_coef_b*((f[idx_f_BB]  + f[idx_f_F]  )*(g[idx_g_BB]  + g[idx_g_F]  )*(h[idx_h_BB]  + h[idx_h_F]  )  +
                                      (f[idx_f_B]   + f[idx_f_FF] )*(g[idx_g_B]   + g[idx_g_FF] )*(h[idx_h_B]   + h[idx_h_FF] )) +
                            d_coef_c*((f[idx_f_BBB] + f[idx_f_F]  )*(g[idx_g_BBB] + g[idx_g_F]  )*(h[idx_h_BBB] + h[idx_h_F]  )  +
                                      (f[idx_f_BB]  + f[idx_f_FF] )*(g[idx_g_BB]  + g[idx_g_FF] )*(h[idx_h_BB]  + h[idx_h_FF] )  +
                                      (f[idx_f_B]   + f[idx_f_FFF])*(g[idx_g_B]   + g[idx_g_FFF])*(h[idx_h_B]   + h[idx_h_FFF])));
                    }
                }
            }
        }
        else if (d_stencil_width == 9)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BBBB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFF  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBBB = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 4 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBB  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_F    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FF   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFF  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFFF = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_BBBB = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 4 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BBB  = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 3 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BB   = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 2 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_B    = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 1 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_F    = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FF   = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 1 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FFF  = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 2 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FFFF = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 3 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_z[idx_face_z] += dt*quarter*(
                            d_coef_a*((f[idx_f_B]    + f[idx_f_F]   )*(g[idx_g_B]    + g[idx_g_F]   )*(h[idx_h_B]    + h[idx_h_F]   )) +
                            d_coef_b*((f[idx_f_BB]   + f[idx_f_F]   )*(g[idx_g_BB]   + g[idx_g_F]   )*(h[idx_h_BB]   + h[idx_h_F]   )  +
                                      (f[idx_f_B]    + f[idx_f_FF]  )*(g[idx_g_B]    + g[idx_g_FF]  )*(h[idx_h_B]    + h[idx_h_FF]  )) +
                            d_coef_c*((f[idx_f_BBB]  + f[idx_f_F]   )*(g[idx_g_BBB]  + g[idx_g_F]   )*(h[idx_h_BBB]  + h[idx_h_F]   )  +
                                      (f[idx_f_BB]   + f[idx_f_FF]  )*(g[idx_g_BB]   + g[idx_g_FF]  )*(h[idx_h_BB]   + h[idx_h_FF]  )  +
                                      (f[idx_f_B]    + f[idx_f_FFF] )*(g[idx_g_B]    + g[idx_g_FFF] )*(h[idx_h_B]    + h[idx_h_FFF] )) +
                            d_coef_d*((f[idx_f_BBBB] + f[idx_f_F]   )*(g[idx_g_BBBB] + g[idx_g_F]   )*(h[idx_h_BBBB] + h[idx_h_F]   )  +
                                      (f[idx_f_BBB]  + f[idx_f_FF]  )*(g[idx_g_BBB]  + g[idx_g_FF]  )*(h[idx_h_BBB]  + h[idx_h_FF]  )  +
                                      (f[idx_f_BB]   + f[idx_f_FFF] )*(g[idx_g_BB]   + g[idx_g_FFF] )*(h[idx_h_BB]   + h[idx_h_FFF] )  +
                                      (f[idx_f_B]    + f[idx_f_FFFF])*(g[idx_g_B]    + g[idx_g_FFFF])*(h[idx_h_B]    + h[idx_h_FFFF])));
                    }
                }
            }
        }
        else if (d_stencil_width == 11)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BBBBB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 5 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBB  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFF   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFF  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFFF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBBBB = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 5 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBBB  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 4 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBB   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B     = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_F     = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FF    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFF   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFFF  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFFFF = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 4 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_BBBBB = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 5 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BBBB  = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 4 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BBB   = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 3 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BB    = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 2 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_B     = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 1 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_F     = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FF    = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 1 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FFF   = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 2 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FFFF  = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 3 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FFFFF = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 4 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_z[idx_face_z] += dt*quarter*(
                            d_coef_a*((f[idx_f_B]     + f[idx_f_F]    )*(g[idx_g_B]     + g[idx_g_F]    )*(h[idx_h_B]     + h[idx_h_F]    )) +
                            d_coef_b*((f[idx_f_BB]    + f[idx_f_F]    )*(g[idx_g_BB]    + g[idx_g_F]    )*(h[idx_h_BB]    + h[idx_h_F]    )  +
                                      (f[idx_f_B]     + f[idx_f_FF]   )*(g[idx_g_B]     + g[idx_g_FF]   )*(h[idx_h_B]     + h[idx_h_FF]   )) +
                            d_coef_c*((f[idx_f_BBB]   + f[idx_f_F]    )*(g[idx_g_BBB]   + g[idx_g_F]    )*(h[idx_h_BBB]   + h[idx_h_F]    )  +
                                      (f[idx_f_BB]    + f[idx_f_FF]   )*(g[idx_g_BB]    + g[idx_g_FF]   )*(h[idx_h_BB]    + h[idx_h_FF]   )  +
                                      (f[idx_f_B]     + f[idx_f_FFF]  )*(g[idx_g_B]     + g[idx_g_FFF]  )*(h[idx_h_B]     + h[idx_h_FFF]  )) +
                            d_coef_d*((f[idx_f_BBBB]  + f[idx_f_F]    )*(g[idx_g_BBBB]  + g[idx_g_F]    )*(h[idx_h_BBBB]  + h[idx_h_F]    )  +
                                      (f[idx_f_BBB]   + f[idx_f_FF]   )*(g[idx_g_BBB]   + g[idx_g_FF]   )*(h[idx_h_BBB]   + h[idx_h_FF]   )  +
                                      (f[idx_f_BB]    + f[idx_f_FFF]  )*(g[idx_g_BB]    + g[idx_g_FFF]  )*(h[idx_h_BB]    + h[idx_h_FFF]  )  +
                                      (f[idx_f_B]     + f[idx_f_FFFF] )*(g[idx_g_B]     + g[idx_g_FFFF] )*(h[idx_h_B]     + h[idx_h_FFFF] )) +
                            d_coef_e*((f[idx_f_BBBBB] + f[idx_f_F]    )*(g[idx_g_BBBBB] + g[idx_g_F]    )*(h[idx_h_BBBBB] + h[idx_h_F]    )  +
                                      (f[idx_f_BBBB]  + f[idx_f_FF]   )*(g[idx_g_BBBB]  + g[idx_g_FF]   )*(h[idx_h_BBBB]  + h[idx_h_FF]   )  +
                                      (f[idx_f_BBB]   + f[idx_f_FFF]  )*(g[idx_g_BBB]   + g[idx_g_FFF]  )*(h[idx_h_BBB]   + h[idx_h_FFF]  )  +
                                      (f[idx_f_BB]    + f[idx_f_FFFF] )*(g[idx_g_BB]    + g[idx_g_FFFF] )*(h[idx_h_BB]    + h[idx_h_FFFF] )  +
                                      (f[idx_f_B]     + f[idx_f_FFFFF])*(g[idx_g_B]     + g[idx_g_FFFFF])*(h[idx_h_B]     + h[idx_h_FFFFF])));
                    }
                }
            }
        }
        else if (d_stencil_width == 13)
        {
            for (int k = 0; k < interior_dim_2 + 1; k++)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx_face_z = i +
                            j*interior_dim_0 + 
                            k*interior_dim_0*interior_dim_1;
                        
                        const int idx_f_BBBBBB = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 6 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBBB  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 5 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBBB   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BBB    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_BB     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_B      = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k - 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_F      = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FF     = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 1 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFF    = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 2 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFF   = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 3 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFFF  = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 4 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_f_FFFFFF = (i + num_ghosts_0_f) +
                            (j + num_ghosts_1_f)*ghostcell_dim_0_f +
                            (k + 5 + num_ghosts_2_f)*ghostcell_dim_0_f*
                                ghostcell_dim_1_f;
                        
                        const int idx_g_BBBBBB = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 6 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBBBB  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 5 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBBB   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 4 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BBB    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_BB     = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_B      = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k - 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_F      = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FF     = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 1 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFF    = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 2 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFFF   = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 3 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFFFF  = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 4 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_g_FFFFFF = (i + num_ghosts_0_g) +
                            (j + num_ghosts_1_g)*ghostcell_dim_0_g +
                            (k + 5 + num_ghosts_2_g)*ghostcell_dim_0_g*
                                ghostcell_dim_1_g;
                        
                        const int idx_h_BBBBBB = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 6 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BBBBB  = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 5 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BBBB   = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 4 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BBB    = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 3 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_BB     = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 2 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_B      = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k - 1 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_F      = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FF     = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 1 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FFF    = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 2 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FFFF   = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 3 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FFFFF  = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 4 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        const int idx_h_FFFFFF = (i + num_ghosts_0_h) +
                            (j + num_ghosts_1_h)*ghostcell_dim_0_h +
                            (k + 5 + num_ghosts_2_h)*ghostcell_dim_0_h*
                                ghostcell_dim_1_h;
                        
                        F_face_z[idx_face_z] += dt*quarter*(
                            d_coef_a*((f[idx_f_B]      + f[idx_f_F]     )*(g[idx_g_B]      + g[idx_g_F]     )*(h[idx_h_B]      + h[idx_h_F]     )) +
                            d_coef_b*((f[idx_f_BB]     + f[idx_f_F]     )*(g[idx_g_BB]     + g[idx_g_F]     )*(h[idx_h_BB]     + h[idx_h_F]     )  +
                                      (f[idx_f_B]      + f[idx_f_FF]    )*(g[idx_g_B]      + g[idx_g_FF]    )*(h[idx_h_B]      + h[idx_h_FF]    )) +
                            d_coef_c*((f[idx_f_BBB]    + f[idx_f_F]     )*(g[idx_g_BBB]    + g[idx_g_F]     )*(h[idx_h_BBB]    + h[idx_h_F]     )  +
                                      (f[idx_f_BB]     + f[idx_f_FF]    )*(g[idx_g_BB]     + g[idx_g_FF]    )*(h[idx_h_BB]     + h[idx_h_FF]    )  +
                                      (f[idx_f_B]      + f[idx_f_FFF]   )*(g[idx_g_B]      + g[idx_g_FFF]   )*(h[idx_h_B]      + h[idx_h_FFF]   )) +
                            d_coef_d*((f[idx_f_BBBB]   + f[idx_f_F]     )*(g[idx_g_BBBB]   + g[idx_g_F]     )*(h[idx_h_BBBB]   + h[idx_h_F]     )  +
                                      (f[idx_f_BBB]    + f[idx_f_FF]    )*(g[idx_g_BBB]    + g[idx_g_FF]    )*(h[idx_h_BBB]    + h[idx_h_FF]    )  +
                                      (f[idx_f_BB]     + f[idx_f_FFF]   )*(g[idx_g_BB]     + g[idx_g_FFF]   )*(h[idx_h_BB]     + h[idx_h_FFF]   )  +
                                      (f[idx_f_B]      + f[idx_f_FFFF]  )*(g[idx_g_B]      + g[idx_g_FFFF]  )*(h[idx_h_B]      + h[idx_h_FFFF]  )) +
                            d_coef_e*((f[idx_f_BBBBB]  + f[idx_f_F]     )*(g[idx_g_BBBBB]  + g[idx_g_F]     )*(h[idx_h_BBBBB]  + h[idx_h_F]     )  +
                                      (f[idx_f_BBBB]   + f[idx_f_FF]    )*(g[idx_g_BBBB]   + g[idx_g_FF]    )*(h[idx_h_BBBB]   + h[idx_h_FF]    )  +
                                      (f[idx_f_BBB]    + f[idx_f_FFF]   )*(g[idx_g_BBB]    + g[idx_g_FFF]   )*(h[idx_h_BBB]    + h[idx_h_FFF]   )  +
                                      (f[idx_f_BB]     + f[idx_f_FFFF]  )*(g[idx_g_BB]     + g[idx_g_FFFF]  )*(h[idx_h_BB]     + h[idx_h_FFFF]  )  +
                                      (f[idx_f_B]      + f[idx_f_FFFFF] )*(g[idx_g_B]      + g[idx_g_FFFFF] )*(h[idx_h_B]      + h[idx_h_FFFFF] )) +
                            d_coef_f*((f[idx_f_BBBBBB] + f[idx_f_F]     )*(g[idx_g_BBBBBB] + g[idx_g_F]     )*(h[idx_h_BBBBBB] + h[idx_h_F]     )  +
                                      (f[idx_f_BBBBB]  + f[idx_f_FF]    )*(g[idx_g_BBBBB]  + g[idx_g_FF]    )*(h[idx_h_BBBBB]  + h[idx_h_FF]    )  +
                                      (f[idx_f_BBBB]   + f[idx_f_FFF]   )*(g[idx_g_BBBB]   + g[idx_g_FFF]   )*(h[idx_h_BBBB]   + h[idx_h_FFF]   )  +
                                      (f[idx_f_BBB]    + f[idx_f_FFFF]  )*(g[idx_g_BBB]    + g[idx_g_FFFF]  )*(h[idx_h_BBB]    + h[idx_h_FFFF]  )  +
                                      (f[idx_f_BB]     + f[idx_f_FFFFF] )*(g[idx_g_BB]     + g[idx_g_FFFFF] )*(h[idx_h_BB]     + h[idx_h_FFFFF] )  +
                                      (f[idx_f_B]      + f[idx_f_FFFFFF])*(g[idx_g_B]      + g[idx_g_FFFFFF])*(h[idx_h_B]      + h[idx_h_FFFFFF])));
                    }
                }
            }
        }
    }
}


/*
 * Add source terms to the advection equations of volume fractions.
 * (for five-equation model by Allaire et al.)
 */
void
ConvectiveFluxReconstructorKEP::addSourceTermsToVolumeFractionEquations(
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_source,
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity,
    HAMERS_SHARED_PTR<pdat::CellData<double> > data_volume_fractions,
    const double* const dx,
    const double dt) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_has_advective_eqn_form);
#endif
    
    const hier::Box interior_box = data_source->getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(data_velocity->getBox().isSpatiallyEqual(interior_box));
    TBOX_ASSERT(data_volume_fractions->getBox().isSpatiallyEqual(interior_box));
#endif
    
    const int num_species = data_volume_fractions->getDepth();
    
    // Get the dimensions of the ghost cell boxes.
    
    const hier::Box ghost_box_velocity = data_velocity->getGhostBox();
    const hier::IntVector ghostcell_dims_velocity = ghost_box_velocity.numberCells();
    const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
    
    const hier::Box ghost_box_volume_fractions = data_volume_fractions->getGhostBox();
    const hier::IntVector ghostcell_dims_volume_fractions = ghost_box_volume_fractions.numberCells();
    const hier::IntVector num_ghosts_volume_fractions = data_volume_fractions->getGhostCellWidth();
    
    /*
     * Get the pointers to the data.
     */
    
    std::vector<double*> Z;
    Z.reserve(num_species);
    for (int si = 0; si < num_species; si++)
    {
        Z.push_back(data_volume_fractions->getPointer(si));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimension.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        /*
         * Get the numbers of ghost cells.
         */
        
        const int num_ghosts_0_velocity         = num_ghosts_velocity[0];
        const int num_ghosts_0_volume_fractions = num_ghosts_volume_fractions[0];
        
        /*
         * Get the pointers to velocity component.
         */
        
        double* u = data_velocity->getPointer(0);
        
        for (int si = 0; si < num_species - 1; si++)
        {
            const int ei = num_species + d_dim.getValue() + 1 + si;
            
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_eqn_form[ei] == EQN_FORM::ADVECTIVE);
#endif
            
            double* S = data_source->getPointer(ei);
            
            if (d_stencil_width == 3)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    
                    const int idx_cell_source = i;
                    
                    const int idx_cell_velocity = i + num_ghosts_0_velocity;
                    
                    const int idx_cell_volume_fractions_x_L = i - 1 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_R = i + 1 + num_ghosts_0_volume_fractions;
                    
                    S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                        d_coef_a*(Z[si][idx_cell_volume_fractions_x_R] - Z[si][idx_cell_volume_fractions_x_L])
                        )/dx[0]
                        );
                }
            }
            else if (d_stencil_width == 5)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    
                    const int idx_cell_source = i;
                    
                    const int idx_cell_velocity = i + num_ghosts_0_velocity;
                    
                    const int idx_cell_volume_fractions_x_LL = i - 2 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_L  = i - 1 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_R  = i + 1 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RR = i + 2 + num_ghosts_0_volume_fractions;
                    
                    S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                        d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]  - Z[si][idx_cell_volume_fractions_x_L]) +
                        d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR] - Z[si][idx_cell_volume_fractions_x_LL])
                        )/dx[0]
                        );
                }
            }
            else if (d_stencil_width == 7)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    
                    const int idx_cell_source = i;
                    
                    const int idx_cell_velocity = i + num_ghosts_0_velocity;
                    
                    const int idx_cell_volume_fractions_x_LLL = i - 3 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_LL  = i - 2 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_L   = i - 1 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_R   = i + 1 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RR  = i + 2 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RRR = i + 3 + num_ghosts_0_volume_fractions;
                    
                    S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                        d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]   - Z[si][idx_cell_volume_fractions_x_L]) +
                        d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR]  - Z[si][idx_cell_volume_fractions_x_LL]) +
                        d_coef_c*(Z[si][idx_cell_volume_fractions_x_RRR] - Z[si][idx_cell_volume_fractions_x_LLL])
                        )/dx[0]
                        );
                }
            }
            else if (d_stencil_width == 9)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    
                    const int idx_cell_source = i;
                    
                    const int idx_cell_velocity = i + num_ghosts_0_velocity;
                    
                    const int idx_cell_volume_fractions_x_LLLL = i - 4 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_LLL  = i - 3 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_LL   = i - 2 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_L    = i - 1 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_R    = i + 1 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RR   = i + 2 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RRR  = i + 3 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RRRR = i + 4 + num_ghosts_0_volume_fractions;
                    
                    S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                        d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]    - Z[si][idx_cell_volume_fractions_x_L]) +
                        d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR]   - Z[si][idx_cell_volume_fractions_x_LL]) +
                        d_coef_c*(Z[si][idx_cell_volume_fractions_x_RRR]  - Z[si][idx_cell_volume_fractions_x_LLL]) +
                        d_coef_d*(Z[si][idx_cell_volume_fractions_x_RRRR] - Z[si][idx_cell_volume_fractions_x_LLLL])
                        )/dx[0]
                        );
                }
            }
            else if (d_stencil_width == 11)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    
                    const int idx_cell_source = i;
                    
                    const int idx_cell_velocity = i + num_ghosts_0_velocity;
                    
                    const int idx_cell_volume_fractions_x_LLLLL = i - 5 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_LLLL  = i - 4 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_LLL   = i - 3 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_LL    = i - 2 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_L     = i - 1 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_R     = i + 1 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RR    = i + 2 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RRR   = i + 3 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RRRR  = i + 4 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RRRRR = i + 5 + num_ghosts_0_volume_fractions;
                    
                    S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                        d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]     - Z[si][idx_cell_volume_fractions_x_L]) +
                        d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR]    - Z[si][idx_cell_volume_fractions_x_LL]) +
                        d_coef_c*(Z[si][idx_cell_volume_fractions_x_RRR]   - Z[si][idx_cell_volume_fractions_x_LLL]) +
                        d_coef_d*(Z[si][idx_cell_volume_fractions_x_RRRR]  - Z[si][idx_cell_volume_fractions_x_LLLL]) +
                        d_coef_e*(Z[si][idx_cell_volume_fractions_x_RRRRR] - Z[si][idx_cell_volume_fractions_x_LLLLL])
                        )/dx[0]
                        );
                }
            }
            else if (d_stencil_width == 13)
            {
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    
                    const int idx_cell_source = i;
                    
                    const int idx_cell_velocity = i + num_ghosts_0_velocity;
                    
                    const int idx_cell_volume_fractions_x_LLLLLL = i - 6 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_LLLLL  = i - 5 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_LLLL   = i - 4 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_LLL    = i - 3 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_LL     = i - 2 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_L      = i - 1 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_R      = i + 1 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RR     = i + 2 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RRR    = i + 3 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RRRR   = i + 4 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RRRRR  = i + 5 + num_ghosts_0_volume_fractions;
                    const int idx_cell_volume_fractions_x_RRRRRR = i + 6 + num_ghosts_0_volume_fractions;
                    
                    S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                        d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]      - Z[si][idx_cell_volume_fractions_x_L]) +
                        d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR]     - Z[si][idx_cell_volume_fractions_x_LL]) +
                        d_coef_c*(Z[si][idx_cell_volume_fractions_x_RRR]    - Z[si][idx_cell_volume_fractions_x_LLL]) +
                        d_coef_d*(Z[si][idx_cell_volume_fractions_x_RRRR]   - Z[si][idx_cell_volume_fractions_x_LLLL]) +
                        d_coef_e*(Z[si][idx_cell_volume_fractions_x_RRRRR]  - Z[si][idx_cell_volume_fractions_x_LLLLL]) +
                        d_coef_f*(Z[si][idx_cell_volume_fractions_x_RRRRRR] - Z[si][idx_cell_volume_fractions_x_LLLLLL])
                        )/dx[0]
                        );
                }
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
        
        const int num_ghosts_0_velocity = num_ghosts_velocity[0];
        const int num_ghosts_1_velocity = num_ghosts_velocity[1];
        const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
        
        const int num_ghosts_0_volume_fractions = num_ghosts_volume_fractions[0];
        const int num_ghosts_1_volume_fractions = num_ghosts_volume_fractions[1];
        const int ghostcell_dim_0_volume_fractions = ghostcell_dims_volume_fractions[0];
        
        /*
         * Get the pointers to velocity components.
         */
        
        double* u = data_velocity->getPointer(0);
        double* v = data_velocity->getPointer(1);
        
        for (int si = 0; si < num_species - 1; si++)
        {
            const int ei = num_species + d_dim.getValue() + 1 + si;
            
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_eqn_form[ei] == EQN_FORM::ADVECTIVE);
#endif
            
            double* S = data_source->getPointer(ei);
            
            if (d_stencil_width == 3)
            {
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        
                        const int idx_cell_source = i +
                            j*interior_dim_0;
                        
                        const int idx_cell_velocity = (i + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_cell_volume_fractions_x_L = (i - 1 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_R = (i + 1 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_B = (i + num_ghosts_0_volume_fractions) +
                            (j - 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_T = (i + num_ghosts_0_volume_fractions) +
                            (j + 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                            d_coef_a*(Z[si][idx_cell_volume_fractions_x_R] - Z[si][idx_cell_volume_fractions_x_L])
                            )/dx[0] +
                            dt*v[idx_cell_velocity]*(
                            d_coef_a*(Z[si][idx_cell_volume_fractions_y_T] - Z[si][idx_cell_volume_fractions_y_B])
                            )/dx[1]
                            );
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
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        
                        const int idx_cell_source = i +
                            j*interior_dim_0;
                        
                        const int idx_cell_velocity = (i + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_cell_volume_fractions_x_LL = (i - 2 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_L  = (i - 1 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_R  = (i + 1 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RR = (i + 2 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BB = (i + num_ghosts_0_volume_fractions) +
                            (j - 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_B  = (i + num_ghosts_0_volume_fractions) +
                            (j - 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_T  = (i + num_ghosts_0_volume_fractions) +
                            (j + 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TT = (i + num_ghosts_0_volume_fractions) +
                            (j + 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                            d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]  - Z[si][idx_cell_volume_fractions_x_L]) +
                            d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR] - Z[si][idx_cell_volume_fractions_x_LL])
                            )/dx[0] +
                            dt*v[idx_cell_velocity]*(
                            d_coef_a*(Z[si][idx_cell_volume_fractions_y_T]  - Z[si][idx_cell_volume_fractions_y_B]) +
                            d_coef_b*(Z[si][idx_cell_volume_fractions_y_TT] - Z[si][idx_cell_volume_fractions_y_BB])
                            )/dx[1]
                            );
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
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        
                        const int idx_cell_source = i +
                            j*interior_dim_0;
                        
                        const int idx_cell_velocity = (i + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_cell_volume_fractions_x_LLL = (i - 3 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_LL  = (i - 2 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_L   = (i - 1 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_R   = (i + 1 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RR  = (i + 2 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RRR = (i + 3 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BBB = (i + num_ghosts_0_volume_fractions) +
                            (j - 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BB  = (i + num_ghosts_0_volume_fractions) +
                            (j - 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_B   = (i + num_ghosts_0_volume_fractions) +
                            (j - 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_T   = (i + num_ghosts_0_volume_fractions) +
                            (j + 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TT  = (i + num_ghosts_0_volume_fractions) +
                            (j + 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TTT = (i + num_ghosts_0_volume_fractions) +
                            (j + 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                            d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]   - Z[si][idx_cell_volume_fractions_x_L]) +
                            d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR]  - Z[si][idx_cell_volume_fractions_x_LL]) +
                            d_coef_c*(Z[si][idx_cell_volume_fractions_x_RRR] - Z[si][idx_cell_volume_fractions_x_LLL])
                            )/dx[0] +
                            dt*v[idx_cell_velocity]*(
                            d_coef_a*(Z[si][idx_cell_volume_fractions_y_T]   - Z[si][idx_cell_volume_fractions_y_B]) +
                            d_coef_b*(Z[si][idx_cell_volume_fractions_y_TT]  - Z[si][idx_cell_volume_fractions_y_BB]) +
                            d_coef_c*(Z[si][idx_cell_volume_fractions_y_TTT] - Z[si][idx_cell_volume_fractions_y_BBB])
                            )/dx[1]
                            );
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
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        
                        const int idx_cell_source = i +
                            j*interior_dim_0;
                        
                        const int idx_cell_velocity = (i + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_cell_volume_fractions_x_LLLL = (i - 4 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_LLL  = (i - 3 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_LL   = (i - 2 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_L    = (i - 1 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_R    = (i + 1 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RR   = (i + 2 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RRR  = (i + 3 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RRRR = (i + 4 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BBBB = (i + num_ghosts_0_volume_fractions) +
                            (j - 4 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BBB  = (i + num_ghosts_0_volume_fractions) +
                            (j - 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BB   = (i + num_ghosts_0_volume_fractions) +
                            (j - 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_B    = (i + num_ghosts_0_volume_fractions) +
                            (j - 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_T    = (i + num_ghosts_0_volume_fractions) +
                            (j + 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TT   = (i + num_ghosts_0_volume_fractions) +
                            (j + 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TTT  = (i + num_ghosts_0_volume_fractions) +
                            (j + 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TTTT = (i + num_ghosts_0_volume_fractions) +
                            (j + 4 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                            d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]    - Z[si][idx_cell_volume_fractions_x_L]) +
                            d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR]   - Z[si][idx_cell_volume_fractions_x_LL]) +
                            d_coef_c*(Z[si][idx_cell_volume_fractions_x_RRR]  - Z[si][idx_cell_volume_fractions_x_LLL]) +
                            d_coef_d*(Z[si][idx_cell_volume_fractions_x_RRRR] - Z[si][idx_cell_volume_fractions_x_LLLL])
                            )/dx[0] +
                            dt*v[idx_cell_velocity]*(
                            d_coef_a*(Z[si][idx_cell_volume_fractions_y_T]    - Z[si][idx_cell_volume_fractions_y_B]) +
                            d_coef_b*(Z[si][idx_cell_volume_fractions_y_TT]   - Z[si][idx_cell_volume_fractions_y_BB]) +
                            d_coef_c*(Z[si][idx_cell_volume_fractions_y_TTT]  - Z[si][idx_cell_volume_fractions_y_BBB]) +
                            d_coef_d*(Z[si][idx_cell_volume_fractions_y_TTTT] - Z[si][idx_cell_volume_fractions_y_BBBB])
                            )/dx[1]
                            );
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
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        
                        const int idx_cell_source = i +
                            j*interior_dim_0;
                        
                        const int idx_cell_velocity = (i + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_cell_volume_fractions_x_LLLLL = (i - 5 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_LLLL  = (i - 4 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_LLL   = (i - 3 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_LL    = (i - 2 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_L     = (i - 1 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_R     = (i + 1 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RR    = (i + 2 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RRR   = (i + 3 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RRRR  = (i + 4 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RRRRR = (i + 5 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BBBBB = (i + num_ghosts_0_volume_fractions) +
                            (j - 5 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BBBB  = (i + num_ghosts_0_volume_fractions) +
                            (j - 4 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BBB   = (i + num_ghosts_0_volume_fractions) +
                            (j - 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BB    = (i + num_ghosts_0_volume_fractions) +
                            (j - 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_B     = (i + num_ghosts_0_volume_fractions) +
                            (j - 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_T     = (i + num_ghosts_0_volume_fractions) +
                            (j + 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TT    = (i + num_ghosts_0_volume_fractions) +
                            (j + 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TTT   = (i + num_ghosts_0_volume_fractions) +
                            (j + 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TTTT  = (i + num_ghosts_0_volume_fractions) +
                            (j + 4 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TTTTT = (i + num_ghosts_0_volume_fractions) +
                            (j + 5 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                            d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]     - Z[si][idx_cell_volume_fractions_x_L]) +
                            d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR]    - Z[si][idx_cell_volume_fractions_x_LL]) +
                            d_coef_c*(Z[si][idx_cell_volume_fractions_x_RRR]   - Z[si][idx_cell_volume_fractions_x_LLL]) +
                            d_coef_d*(Z[si][idx_cell_volume_fractions_x_RRRR]  - Z[si][idx_cell_volume_fractions_x_LLLL]) +
                            d_coef_e*(Z[si][idx_cell_volume_fractions_x_RRRRR] - Z[si][idx_cell_volume_fractions_x_LLLLL])
                            )/dx[0] +
                            dt*v[idx_cell_velocity]*(
                            d_coef_a*(Z[si][idx_cell_volume_fractions_y_T]     - Z[si][idx_cell_volume_fractions_y_B]) +
                            d_coef_b*(Z[si][idx_cell_volume_fractions_y_TT]    - Z[si][idx_cell_volume_fractions_y_BB]) +
                            d_coef_c*(Z[si][idx_cell_volume_fractions_y_TTT]   - Z[si][idx_cell_volume_fractions_y_BBB]) +
                            d_coef_d*(Z[si][idx_cell_volume_fractions_y_TTTT]  - Z[si][idx_cell_volume_fractions_y_BBBB]) +
                            d_coef_e*(Z[si][idx_cell_volume_fractions_y_TTTTT] - Z[si][idx_cell_volume_fractions_y_BBBBB])
                            )/dx[1]
                            );
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
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        
                        const int idx_cell_source = i +
                            j*interior_dim_0;
                        
                        const int idx_cell_velocity = (i + num_ghosts_0_velocity) +
                            (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                        
                        const int idx_cell_volume_fractions_x_LLLLLL = (i - 6 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_LLLLL  = (i - 5 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_LLLL   = (i - 4 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_LLL    = (i - 3 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_LL     = (i - 2 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_L      = (i - 1 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_R      = (i + 1 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RR     = (i + 2 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RRR    = (i + 3 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RRRR   = (i + 4 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RRRRR  = (i + 5 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_x_RRRRRR = (i + 6 + num_ghosts_0_volume_fractions) +
                            (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BBBBBB = (i + num_ghosts_0_volume_fractions) +
                            (j - 6 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BBBBB  = (i + num_ghosts_0_volume_fractions) +
                            (j - 5 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BBBB   = (i + num_ghosts_0_volume_fractions) +
                            (j - 4 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BBB    = (i + num_ghosts_0_volume_fractions) +
                            (j - 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_BB     = (i + num_ghosts_0_volume_fractions) +
                            (j - 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_B      = (i + num_ghosts_0_volume_fractions) +
                            (j - 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_T      = (i + num_ghosts_0_volume_fractions) +
                            (j + 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TT     = (i + num_ghosts_0_volume_fractions) +
                            (j + 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TTT    = (i + num_ghosts_0_volume_fractions) +
                            (j + 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TTTT   = (i + num_ghosts_0_volume_fractions) +
                            (j + 4 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TTTTT  = (i + num_ghosts_0_volume_fractions) +
                            (j + 5 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        const int idx_cell_volume_fractions_y_TTTTTT = (i + num_ghosts_0_volume_fractions) +
                            (j + 6 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions;
                        
                        S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                            d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]      - Z[si][idx_cell_volume_fractions_x_L]) +
                            d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR]     - Z[si][idx_cell_volume_fractions_x_LL]) +
                            d_coef_c*(Z[si][idx_cell_volume_fractions_x_RRR]    - Z[si][idx_cell_volume_fractions_x_LLL]) +
                            d_coef_d*(Z[si][idx_cell_volume_fractions_x_RRRR]   - Z[si][idx_cell_volume_fractions_x_LLLL]) +
                            d_coef_e*(Z[si][idx_cell_volume_fractions_x_RRRRR]  - Z[si][idx_cell_volume_fractions_x_LLLLL]) +
                            d_coef_f*(Z[si][idx_cell_volume_fractions_x_RRRRRR] - Z[si][idx_cell_volume_fractions_x_LLLLLL])
                            )/dx[0] +
                            dt*v[idx_cell_velocity]*(
                            d_coef_a*(Z[si][idx_cell_volume_fractions_y_T]      - Z[si][idx_cell_volume_fractions_y_B]) +
                            d_coef_b*(Z[si][idx_cell_volume_fractions_y_TT]     - Z[si][idx_cell_volume_fractions_y_BB]) +
                            d_coef_c*(Z[si][idx_cell_volume_fractions_y_TTT]    - Z[si][idx_cell_volume_fractions_y_BBB]) +
                            d_coef_d*(Z[si][idx_cell_volume_fractions_y_TTTT]   - Z[si][idx_cell_volume_fractions_y_BBBB]) +
                            d_coef_e*(Z[si][idx_cell_volume_fractions_y_TTTTT]  - Z[si][idx_cell_volume_fractions_y_BBBBB]) +
                            d_coef_f*(Z[si][idx_cell_volume_fractions_y_TTTTTT] - Z[si][idx_cell_volume_fractions_y_BBBBBB])
                            )/dx[1]
                            );
                    }
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
        
        const int num_ghosts_0_velocity = num_ghosts_velocity[0];
        const int num_ghosts_1_velocity = num_ghosts_velocity[1];
        const int num_ghosts_2_velocity = num_ghosts_velocity[2];
        const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
        const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
        
        const int num_ghosts_0_volume_fractions = num_ghosts_volume_fractions[0];
        const int num_ghosts_1_volume_fractions = num_ghosts_volume_fractions[1];
        const int num_ghosts_2_volume_fractions = num_ghosts_volume_fractions[2];
        const int ghostcell_dim_0_volume_fractions = ghostcell_dims_volume_fractions[0];
        const int ghostcell_dim_1_volume_fractions = ghostcell_dims_volume_fractions[1];
        
        /*
         * Get the pointers to velocity components.
         */
        
        double* u = data_velocity->getPointer(0);
        double* v = data_velocity->getPointer(1);
        double* w = data_velocity->getPointer(2);
        
        for (int si = 0; si < num_species - 1; si++)
        {
            const int ei = num_species + d_dim.getValue() + 1 + si;
            
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_eqn_form[ei] == EQN_FORM::ADVECTIVE);
#endif
            
            double* S = data_source->getPointer(ei);
            
            if (d_stencil_width == 3)
            {
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
                            
                            const int idx_cell_source = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_cell_velocity = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_cell_volume_fractions_x_L = (i - 1 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_R = (i + 1 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_B = (i + num_ghosts_0_volume_fractions) +
                                (j - 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_T = (i + num_ghosts_0_volume_fractions) +
                                (j + 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_B = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 1 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_F = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 1 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_x_R] - Z[si][idx_cell_volume_fractions_x_L])
                                )/dx[0] +
                                dt*v[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_y_T] - Z[si][idx_cell_volume_fractions_y_B])
                                )/dx[1] +
                                dt*w[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_z_F] - Z[si][idx_cell_volume_fractions_z_B])
                                )/dx[2]
                                );
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
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            
                            const int idx_cell_source = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_cell_velocity = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_cell_volume_fractions_x_LL = (i - 2 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_L  = (i - 1 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_R  = (i + 1 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RR = (i + 2 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BB = (i + num_ghosts_0_volume_fractions) +
                                (j - 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_B  = (i + num_ghosts_0_volume_fractions) +
                                (j - 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_T  = (i + num_ghosts_0_volume_fractions) +
                                (j + 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TT = (i + num_ghosts_0_volume_fractions) +
                                (j + 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BB = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 2 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_B  = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 1 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_F  = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 1 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FF = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 2 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]  - Z[si][idx_cell_volume_fractions_x_L]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR] - Z[si][idx_cell_volume_fractions_x_LL])
                                )/dx[0] +
                                dt*v[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_y_T]  - Z[si][idx_cell_volume_fractions_y_B]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_y_TT] - Z[si][idx_cell_volume_fractions_y_BB])
                                )/dx[1] +
                                dt*w[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_z_F]  - Z[si][idx_cell_volume_fractions_z_B]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_z_FF] - Z[si][idx_cell_volume_fractions_z_BB])
                                )/dx[2]
                                );
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
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            
                            const int idx_cell_source = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_cell_velocity = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_cell_volume_fractions_x_LLL = (i - 3 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_LL  = (i - 2 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_L   = (i - 1 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_R   = (i + 1 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RR  = (i + 2 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RRR = (i + 3 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BBB = (i + num_ghosts_0_volume_fractions) +
                                (j - 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BB  = (i + num_ghosts_0_volume_fractions) +
                                (j - 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_B   = (i + num_ghosts_0_volume_fractions) +
                                (j - 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_T   = (i + num_ghosts_0_volume_fractions) +
                                (j + 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TT  = (i + num_ghosts_0_volume_fractions) +
                                (j + 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TTT = (i + num_ghosts_0_volume_fractions) +
                                (j + 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BBB = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 3 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BB  = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 2 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_B   = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 1 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_F   = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 1 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FF  = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 2 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FFF = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 3 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]   - Z[si][idx_cell_volume_fractions_x_L]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR]  - Z[si][idx_cell_volume_fractions_x_LL]) +
                                d_coef_c*(Z[si][idx_cell_volume_fractions_x_RRR] - Z[si][idx_cell_volume_fractions_x_LLL])
                                )/dx[0] +
                                dt*v[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_y_T]   - Z[si][idx_cell_volume_fractions_y_B]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_y_TT]  - Z[si][idx_cell_volume_fractions_y_BB]) +
                                d_coef_c*(Z[si][idx_cell_volume_fractions_y_TTT] - Z[si][idx_cell_volume_fractions_y_BBB])
                                )/dx[1] +
                                dt*w[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_z_F]   - Z[si][idx_cell_volume_fractions_z_B]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_z_FF]  - Z[si][idx_cell_volume_fractions_z_BB]) +
                                d_coef_c*(Z[si][idx_cell_volume_fractions_z_FFF] - Z[si][idx_cell_volume_fractions_z_BBB])
                                )/dx[2]
                                );
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
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            
                            const int idx_cell_source = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_cell_velocity = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_cell_volume_fractions_x_LLLL = (i - 4 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_LLL  = (i - 3 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_LL   = (i - 2 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_L    = (i - 1 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_R    = (i + 1 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RR   = (i + 2 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RRR  = (i + 3 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RRRR = (i + 4 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BBBB = (i + num_ghosts_0_volume_fractions) +
                                (j - 4 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BBB  = (i + num_ghosts_0_volume_fractions) +
                                (j - 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BB   = (i + num_ghosts_0_volume_fractions) +
                                (j - 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_B    = (i + num_ghosts_0_volume_fractions) +
                                (j - 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_T    = (i + num_ghosts_0_volume_fractions) +
                                (j + 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TT   = (i + num_ghosts_0_volume_fractions) +
                                (j + 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TTT  = (i + num_ghosts_0_volume_fractions) +
                                (j + 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TTTT = (i + num_ghosts_0_volume_fractions) +
                                (j + 4 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BBBB = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 4 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BBB  = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 3 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BB   = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 2 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_B    = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 1 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_F    = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 1 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FF   = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 2 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FFF  = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 3 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FFFF = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 4 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]    - Z[si][idx_cell_volume_fractions_x_L]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR]   - Z[si][idx_cell_volume_fractions_x_LL]) +
                                d_coef_c*(Z[si][idx_cell_volume_fractions_x_RRR]  - Z[si][idx_cell_volume_fractions_x_LLL]) +
                                d_coef_d*(Z[si][idx_cell_volume_fractions_x_RRRR] - Z[si][idx_cell_volume_fractions_x_LLLL])
                                )/dx[0] +
                                dt*v[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_y_T]    - Z[si][idx_cell_volume_fractions_y_B]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_y_TT]   - Z[si][idx_cell_volume_fractions_y_BB]) +
                                d_coef_c*(Z[si][idx_cell_volume_fractions_y_TTT]  - Z[si][idx_cell_volume_fractions_y_BBB]) +
                                d_coef_d*(Z[si][idx_cell_volume_fractions_y_TTTT] - Z[si][idx_cell_volume_fractions_y_BBBB])
                                )/dx[1] +
                                dt*w[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_z_F]    - Z[si][idx_cell_volume_fractions_z_B]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_z_FF]   - Z[si][idx_cell_volume_fractions_z_BB]) +
                                d_coef_c*(Z[si][idx_cell_volume_fractions_z_FFF]  - Z[si][idx_cell_volume_fractions_z_BBB]) +
                                d_coef_d*(Z[si][idx_cell_volume_fractions_z_FFFF] - Z[si][idx_cell_volume_fractions_z_BBBB])
                                )/dx[2]
                                );
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
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            
                            const int idx_cell_source = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_cell_velocity = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_cell_volume_fractions_x_LLLLL = (i - 5 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_LLLL  = (i - 4 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_LLL   = (i - 3 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_LL    = (i - 2 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_L     = (i - 1 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_R     = (i + 1 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RR    = (i + 2 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RRR   = (i + 3 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RRRR  = (i + 4 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RRRRR = (i + 5 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BBBBB = (i + num_ghosts_0_volume_fractions) +
                                (j - 5 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BBBB  = (i + num_ghosts_0_volume_fractions) +
                                (j - 4 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BBB   = (i + num_ghosts_0_volume_fractions) +
                                (j - 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BB    = (i + num_ghosts_0_volume_fractions) +
                                (j - 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_B     = (i + num_ghosts_0_volume_fractions) +
                                (j - 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_T     = (i + num_ghosts_0_volume_fractions) +
                                (j + 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TT    = (i + num_ghosts_0_volume_fractions) +
                                (j + 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TTT   = (i + num_ghosts_0_volume_fractions) +
                                (j + 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TTTT  = (i + num_ghosts_0_volume_fractions) +
                                (j + 4 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TTTTT = (i + num_ghosts_0_volume_fractions) +
                                (j + 5 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BBBBB = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 5 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BBBB  = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 4 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BBB   = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 3 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BB    = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 2 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_B     = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 1 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_F     = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 1 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FF    = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 2 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FFF   = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 3 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FFFF  = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 4 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FFFFF = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 5 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]     - Z[si][idx_cell_volume_fractions_x_L]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR]    - Z[si][idx_cell_volume_fractions_x_LL]) +
                                d_coef_c*(Z[si][idx_cell_volume_fractions_x_RRR]   - Z[si][idx_cell_volume_fractions_x_LLL]) +
                                d_coef_d*(Z[si][idx_cell_volume_fractions_x_RRRR]  - Z[si][idx_cell_volume_fractions_x_LLLL]) +
                                d_coef_e*(Z[si][idx_cell_volume_fractions_x_RRRRR] - Z[si][idx_cell_volume_fractions_x_LLLLL])
                                )/dx[0] +
                                dt*v[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_y_T]     - Z[si][idx_cell_volume_fractions_y_B]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_y_TT]    - Z[si][idx_cell_volume_fractions_y_BB]) +
                                d_coef_c*(Z[si][idx_cell_volume_fractions_y_TTT]   - Z[si][idx_cell_volume_fractions_y_BBB]) +
                                d_coef_d*(Z[si][idx_cell_volume_fractions_y_TTTT]  - Z[si][idx_cell_volume_fractions_y_BBBB]) +
                                d_coef_e*(Z[si][idx_cell_volume_fractions_y_TTTTT] - Z[si][idx_cell_volume_fractions_y_BBBBB])
                                )/dx[1] +
                                dt*w[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_z_F]     - Z[si][idx_cell_volume_fractions_z_B]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_z_FF]    - Z[si][idx_cell_volume_fractions_z_BB]) +
                                d_coef_c*(Z[si][idx_cell_volume_fractions_z_FFF]   - Z[si][idx_cell_volume_fractions_z_BBB]) +
                                d_coef_d*(Z[si][idx_cell_volume_fractions_z_FFFF]  - Z[si][idx_cell_volume_fractions_z_BBBB]) +
                                d_coef_e*(Z[si][idx_cell_volume_fractions_z_FFFFF] - Z[si][idx_cell_volume_fractions_z_BBBBB])
                                )/dx[2]
                                );
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
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            
                            const int idx_cell_source = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*
                                    interior_dim_1;
                            
                            const int idx_cell_velocity = (i + num_ghosts_0_velocity) +
                                (j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                (k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                    ghostcell_dim_1_velocity;
                            
                            const int idx_cell_volume_fractions_x_LLLLLL = (i - 6 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_LLLLL  = (i - 5 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_LLLL   = (i - 4 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_LLL    = (i - 3 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_LL     = (i - 2 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_L      = (i - 1 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_R      = (i + 1 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RR     = (i + 2 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RRR    = (i + 3 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RRRR   = (i + 4 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RRRRR  = (i + 5 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_x_RRRRRR = (i + 6 + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BBBBBB = (i + num_ghosts_0_volume_fractions) +
                                (j - 6 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BBBBB  = (i + num_ghosts_0_volume_fractions) +
                                (j - 5 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BBBB   = (i + num_ghosts_0_volume_fractions) +
                                (j - 4 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BBB    = (i + num_ghosts_0_volume_fractions) +
                                (j - 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_BB     = (i + num_ghosts_0_volume_fractions) +
                                (j - 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_B      = (i + num_ghosts_0_volume_fractions) +
                                (j - 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_T      = (i + num_ghosts_0_volume_fractions) +
                                (j + 1 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TT     = (i + num_ghosts_0_volume_fractions) +
                                (j + 2 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TTT    = (i + num_ghosts_0_volume_fractions) +
                                (j + 3 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TTTT   = (i + num_ghosts_0_volume_fractions) +
                                (j + 4 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TTTTT  = (i + num_ghosts_0_volume_fractions) +
                                (j + 5 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_y_TTTTTT = (i + num_ghosts_0_volume_fractions) +
                                (j + 6 + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BBBBBB = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 6 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BBBBB  = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 5 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BBBB   = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 4 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BBB    = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 3 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_BB     = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 2 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_B      = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k - 1 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_F      = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 1 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FF     = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 2 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FFF    = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 3 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FFFF   = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 4 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FFFFF  = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 5 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            const int idx_cell_volume_fractions_z_FFFFFF = (i + num_ghosts_0_volume_fractions) +
                                (j + num_ghosts_1_volume_fractions)*ghostcell_dim_0_volume_fractions +
                                (k + 6 + num_ghosts_2_volume_fractions)*ghostcell_dim_0_volume_fractions*
                                    ghostcell_dim_1_volume_fractions;
                            
                            S[idx_cell_source] -= (dt*u[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_x_R]      - Z[si][idx_cell_volume_fractions_x_L]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_x_RR]     - Z[si][idx_cell_volume_fractions_x_LL]) +
                                d_coef_c*(Z[si][idx_cell_volume_fractions_x_RRR]    - Z[si][idx_cell_volume_fractions_x_LLL]) +
                                d_coef_d*(Z[si][idx_cell_volume_fractions_x_RRRR]   - Z[si][idx_cell_volume_fractions_x_LLLL]) +
                                d_coef_e*(Z[si][idx_cell_volume_fractions_x_RRRRR]  - Z[si][idx_cell_volume_fractions_x_LLLLL]) +
                                d_coef_f*(Z[si][idx_cell_volume_fractions_x_RRRRRR] - Z[si][idx_cell_volume_fractions_x_LLLLLL])
                                )/dx[0] +
                                dt*v[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_y_T]      - Z[si][idx_cell_volume_fractions_y_B]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_y_TT]     - Z[si][idx_cell_volume_fractions_y_BB]) +
                                d_coef_c*(Z[si][idx_cell_volume_fractions_y_TTT]    - Z[si][idx_cell_volume_fractions_y_BBB]) +
                                d_coef_d*(Z[si][idx_cell_volume_fractions_y_TTTT]   - Z[si][idx_cell_volume_fractions_y_BBBB]) +
                                d_coef_e*(Z[si][idx_cell_volume_fractions_y_TTTTT]  - Z[si][idx_cell_volume_fractions_y_BBBBB]) +
                                d_coef_f*(Z[si][idx_cell_volume_fractions_y_TTTTTT] - Z[si][idx_cell_volume_fractions_y_BBBBBB])
                                )/dx[1] +
                                dt*w[idx_cell_velocity]*(
                                d_coef_a*(Z[si][idx_cell_volume_fractions_z_F]      - Z[si][idx_cell_volume_fractions_z_B]) +
                                d_coef_b*(Z[si][idx_cell_volume_fractions_z_FF]     - Z[si][idx_cell_volume_fractions_z_BB]) +
                                d_coef_c*(Z[si][idx_cell_volume_fractions_z_FFF]    - Z[si][idx_cell_volume_fractions_z_BBB]) +
                                d_coef_d*(Z[si][idx_cell_volume_fractions_z_FFFF]   - Z[si][idx_cell_volume_fractions_z_BBBB]) +
                                d_coef_e*(Z[si][idx_cell_volume_fractions_z_FFFFF]  - Z[si][idx_cell_volume_fractions_z_BBBBB]) +
                                d_coef_f*(Z[si][idx_cell_volume_fractions_z_FFFFFF] - Z[si][idx_cell_volume_fractions_z_BBBBBB])
                                )/dx[2]
                                );
                        }
                    }
                }
            }
        }
    }
}

