#include "flow/convective_flux_reconstructors/first_order/ConvectiveFluxReconstructorFirstOrderLLF.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

ConvectiveFluxReconstructorFirstOrderLLF::ConvectiveFluxReconstructorFirstOrderLLF(
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
    d_num_conv_ghosts = hier::IntVector::getOne(d_dim);
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorFirstOrderLLF::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorFirstOrderLLF object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorFirstOrderLLF: this = "
       << (ConvectiveFluxReconstructorFirstOrderLLF *)this
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
ConvectiveFluxReconstructorFirstOrderLLF::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putString("d_shock_capturing_scheme", "LLF");
}


/*
 * Compute the convective flux and source due to splitting of convective term on a patch.
 */
void
ConvectiveFluxReconstructorFirstOrderLLF::computeConvectiveFluxAndSourceOnPatch(
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
    
    /*
     * Get the forms of equation and check whether there are advection equations in the system of equations.
     */
    
    const std::vector<EQN_FORM::TYPE> eqn_form = d_flow_model->getEquationsForm();
    
    bool has_advection_eqn = false;
    
    if (std::find(eqn_form.begin(), eqn_form.end(), EQN_FORM::ADVECTIVE) != eqn_form.end())
    {
        has_advection_eqn = true;
    }
    
    // Allocate temporary patch data.
    HAMERS_SHARED_PTR<pdat::SideData<Real> > velocity_intercell;
    
    if (has_advection_eqn)
    {
        velocity_intercell.reset(
            new pdat::SideData<Real>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        if (has_advection_eqn)
        {
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("DENSITY", d_num_conv_ghosts));
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRESSURE", d_num_conv_ghosts));
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
        }
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("MAX_WAVE_SPEED_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        
        d_flow_model->registerDerivedVariables(num_subghosts_of_data);
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Get the pointer to the cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > density;
        HAMERS_SHARED_PTR<pdat::CellData<Real> > pressure;
        HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity;
        
        if (has_advection_eqn)
        {
            density = d_flow_model->getCellData("DENSITY");
            pressure = d_flow_model->getCellData("PRESSURE");
            velocity = d_flow_model->getCellData("VELOCITY");
        }
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > max_wave_speed_x = d_flow_model->getCellData("MAX_WAVE_SPEED_X");
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(1);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        
        hier::IntVector num_subghosts_density(d_dim);
        hier::IntVector num_subghosts_pressure(d_dim);
        hier::IntVector num_subghosts_velocity(d_dim);
        
        hier::IntVector subghostcell_dims_density(d_dim);
        hier::IntVector subghostcell_dims_pressure(d_dim);
        hier::IntVector subghostcell_dims_velocity(d_dim);
        
        if (has_advection_eqn)
        {
            num_subghosts_density = density->getGhostCellWidth();
            subghostcell_dims_density = density->getGhostBox().numberCells();
            
            num_subghosts_pressure = pressure->getGhostCellWidth();
            subghostcell_dims_pressure = pressure->getGhostBox().numberCells();
            
            num_subghosts_velocity = velocity->getGhostCellWidth();
            subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
        }
        
        hier::IntVector num_subghosts_max_wave_speed_x = max_wave_speed_x->getGhostCellWidth();
        hier::IntVector subghostcell_dims_max_wave_speed_x = max_wave_speed_x->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
        
        Real* rho = nullptr;
        Real* p = nullptr;
        Real* u = nullptr;
        
        if (has_advection_eqn)
        {
            rho = density->getPointer(0);
            p = pressure->getPointer(0);
            u = velocity->getPointer(0);
        }
        
        Real* max_lambda_x = max_wave_speed_x->getPointer(0);
        
        std::vector<Real*> F_x_node;
        F_x_node.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
        }
        
        /*
         * Get the pointers to the conservative variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
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
        
        for (int i = 0; i < interior_dims[0] + 1; i++)
        {
            // Compute the linear indices.
            const int idx_midpoint_x = i;
            const int idx_L_max_wave_speed_x = i - 1 + num_subghosts_max_wave_speed_x[0];
            const int idx_R_max_wave_speed_x = i + num_subghosts_max_wave_speed_x[0];
            const int idx_L_convective_flux_x = i - 1 + num_subghosts_convective_flux_x[0];
            const int idx_R_convective_flux_x = i + num_subghosts_convective_flux_x[0];
            
            const Real alpha_x = std::max(max_lambda_x[idx_L_max_wave_speed_x], max_lambda_x[idx_R_max_wave_speed_x]);
            
            // Compute the fluxes.
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                // Compute the linear indices.
                const int idx_L_conservative_var = i - 1 + num_subghosts_conservative_var[ei][0];
                const int idx_R_conservative_var = i + num_subghosts_conservative_var[ei][0];
                
                convective_flux->getPointer(0, ei)[idx_midpoint_x] = Real(1)/Real(2)*Real(dt)*(
                    F_x_node[ei][idx_L_convective_flux_x] + F_x_node[ei][idx_R_convective_flux_x] -
                        alpha_x*(Q[ei][idx_R_conservative_var] - Q[ei][idx_L_conservative_var]));
            }
            
            if (has_advection_eqn)
            {
                // Compute the linear indices.
                const int idx_L_density = i - 1 + num_subghosts_density[0];
                const int idx_R_density = i + num_subghosts_density[0];
                const int idx_L_pressure = i - 1 + num_subghosts_pressure[0];
                const int idx_R_pressure = i + num_subghosts_pressure[0];
                const int idx_L_velocity = i - 1 + num_subghosts_velocity[0];
                const int idx_R_velocity = i + num_subghosts_velocity[0];
                
                velocity_intercell->getPointer(0, 0)[idx_midpoint_x] = (p[idx_R_pressure] - p[idx_L_pressure] +
                    rho[idx_L_density]*u[idx_L_velocity]*(-alpha_x - u[idx_L_velocity]) -
                        rho[idx_R_density]*u[idx_R_velocity]*(alpha_x - u[idx_R_velocity]))/
                            (rho[idx_L_density]*(-alpha_x - u[idx_L_velocity]) -
                                rho[idx_R_density]*(alpha_x - u[idx_R_velocity]));
            }
        }
        
        /*
         * Compute the source.
         */
        
        if (has_advection_eqn)
        {
            for (int ei = 0; ei < d_num_eqn; ei ++)
            {
                if (eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute the linear indices. 
                        const int idx_cell_wghost = i + num_subghosts_conservative_var[ei][0];
                        const int idx_cell_nghost = i;
                        const int idx_midpoint_x_L = i;
                        const int idx_midpoint_x_R = i + 1;
                        
                        const Real& u_L = velocity_intercell->getPointer(0, 0)[idx_midpoint_x_L];
                        const Real& u_R = velocity_intercell->getPointer(0, 0)[idx_midpoint_x_R];
                        
                        S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(u_R - u_L)/Real(dx[0]);
                    }
                }
            }
        }
        
        /*
         * Unregister the patch and data of all registered derived cell variables in the flow model.
         */
        
        d_flow_model->unregisterPatch();
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        if (has_advection_eqn)
        {
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("DENSITY", d_num_conv_ghosts));
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRESSURE", d_num_conv_ghosts));
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
        }
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("MAX_WAVE_SPEED_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("MAX_WAVE_SPEED_Y", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
        
        d_flow_model->registerDerivedVariables(num_subghosts_of_data);
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Get the pointer to the cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > density;
        HAMERS_SHARED_PTR<pdat::CellData<Real> > pressure;
        HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity;
        
        if (has_advection_eqn)
        {
            density = d_flow_model->getCellData("DENSITY");
            pressure = d_flow_model->getCellData("PRESSURE");
            velocity = d_flow_model->getCellData("VELOCITY");
        }
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > max_wave_speed_x = d_flow_model->getCellData("MAX_WAVE_SPEED_X");
        HAMERS_SHARED_PTR<pdat::CellData<Real> > max_wave_speed_y = d_flow_model->getCellData("MAX_WAVE_SPEED_Y");
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(2);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
        
        hier::IntVector num_subghosts_density(d_dim);
        hier::IntVector num_subghosts_pressure(d_dim);
        hier::IntVector num_subghosts_velocity(d_dim);
        
        hier::IntVector subghostcell_dims_density(d_dim);
        hier::IntVector subghostcell_dims_pressure(d_dim);
        hier::IntVector subghostcell_dims_velocity(d_dim);
        
        if (has_advection_eqn)
        {
            num_subghosts_density = density->getGhostCellWidth();
            subghostcell_dims_density = density->getGhostBox().numberCells();
            
            num_subghosts_pressure = pressure->getGhostCellWidth();
            subghostcell_dims_pressure = pressure->getGhostBox().numberCells();
            
            num_subghosts_velocity = velocity->getGhostCellWidth();
            subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
        }
        
        hier::IntVector num_subghosts_max_wave_speed_x = max_wave_speed_x->getGhostCellWidth();
        hier::IntVector subghostcell_dims_max_wave_speed_x = max_wave_speed_x->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_max_wave_speed_y = max_wave_speed_y->getGhostCellWidth();
        hier::IntVector subghostcell_dims_max_wave_speed_y = max_wave_speed_y->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
        
        Real* rho = nullptr;
        Real* p = nullptr;
        Real* u = nullptr;
        Real* v = nullptr;
        
        if (has_advection_eqn)
        {
            rho = density->getPointer(0);
            p = pressure->getPointer(0);
            u = velocity->getPointer(0);
            v = velocity->getPointer(1);
        }
        
        Real* max_lambda_x = max_wave_speed_x->getPointer(0);
        Real* max_lambda_y = max_wave_speed_y->getPointer(0);
        
        std::vector<Real*> F_x_node;
        std::vector<Real*> F_y_node;
        F_x_node.reserve(d_num_eqn);
        F_y_node.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
            F_y_node.push_back(convective_flux_node[1]->getPointer(ei));
        }
        
        /*
         * Get the pointers to the conservative variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
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

        for (int j = 0; j < interior_dims[1]; j++)
        {
            for (int i = 0; i < interior_dims[0] + 1; i++)
            {
                // Compute the linear indices.
                const int idx_midpoint_x = i +
                    j*(interior_dims[0] + 1);
                
                const int idx_L_max_wave_speed_x = (i - 1 + num_subghosts_max_wave_speed_x[0]) +
                    (j + num_subghosts_max_wave_speed_x[1])*subghostcell_dims_max_wave_speed_x[0];
                
                const int idx_R_max_wave_speed_x = (i + num_subghosts_max_wave_speed_x[0]) +
                    (j + num_subghosts_max_wave_speed_x[1])*subghostcell_dims_max_wave_speed_x[0];
                
                const int idx_L_convective_flux_x = (i - 1 + num_subghosts_convective_flux_x[0]) +
                    (j + num_subghosts_convective_flux_x[1])*subghostcell_dims_convective_flux_x[0];
                
                const int idx_R_convective_flux_x = (i + num_subghosts_convective_flux_x[0]) +
                    (j + num_subghosts_convective_flux_x[1])*subghostcell_dims_convective_flux_x[0];
                
                const Real alpha_x = std::max(max_lambda_x[idx_L_max_wave_speed_x], max_lambda_x[idx_R_max_wave_speed_x]);
                
                // Compute the fluxes.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    // Compute the linear indices.
                    const int idx_L_conservative_var = (i - 1 + num_subghosts_conservative_var[ei][0]) +
                        (j + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0];
                    
                    const int idx_R_conservative_var = (i + num_subghosts_conservative_var[ei][0]) +
                        (j + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0];
                    
                    convective_flux->getPointer(0, ei)[idx_midpoint_x] = Real(1)/Real(2)*Real(dt)*(
                        F_x_node[ei][idx_L_convective_flux_x] + F_x_node[ei][idx_R_convective_flux_x] -
                            alpha_x*(Q[ei][idx_R_conservative_var] - Q[ei][idx_L_conservative_var]));
                }
                
                if (has_advection_eqn)
                {
                    // Compute the linear indices.
                    const int idx_L_density = (i - 1 + num_subghosts_density[0]) +
                        (j + num_subghosts_density[1])*subghostcell_dims_density[0];
                    
                    const int idx_R_density = (i + num_subghosts_density[0]) +
                        (j + num_subghosts_density[1])*subghostcell_dims_density[0];
                    
                    const int idx_L_pressure = (i - 1 + num_subghosts_pressure[0]) +
                        (j + num_subghosts_pressure[1])*subghostcell_dims_pressure[0];
                    
                    const int idx_R_pressure = (i + num_subghosts_pressure[0]) +
                        (j + num_subghosts_pressure[1])*subghostcell_dims_pressure[0];
                    
                    const int idx_L_velocity = (i - 1 + num_subghosts_velocity[0]) +
                        (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                    
                    const int idx_R_velocity = (i + num_subghosts_velocity[0]) +
                        (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                    
                    velocity_intercell->getPointer(0, 0)[idx_midpoint_x] = (p[idx_R_pressure] - p[idx_L_pressure] +
                        rho[idx_L_density]*u[idx_L_velocity]*(-alpha_x - u[idx_L_velocity]) -
                            rho[idx_R_density]*u[idx_R_velocity]*(alpha_x - u[idx_R_velocity]))/
                                (rho[idx_L_density]*(-alpha_x - u[idx_L_velocity]) -
                                    rho[idx_R_density]*(alpha_x - u[idx_R_velocity]));
                    
                    velocity_intercell->getPointer(0, 1)[idx_midpoint_x] =
                        (rho[idx_L_density]*v[idx_L_velocity]*(-alpha_x - u[idx_L_velocity]) -
                            rho[idx_R_density]*v[idx_R_velocity]*(alpha_x - u[idx_R_velocity]))/
                                (rho[idx_L_density]*(-alpha_x - u[idx_L_velocity]) -
                                    rho[idx_R_density]*(alpha_x - u[idx_R_velocity]));
                }
            }
        }
        
        /*
         * Compute the fluxes in the y direction.
         */
        
        for (int j = 0; j < interior_dims[1] + 1; j++)
        {
            for (int i = 0; i < interior_dims[0]; i++)
            {
                // Compute the linear indices.
                const int idx_midpoint_y = i +
                    j*interior_dims[0];
                
                const int idx_B_max_wave_speed_y = (i + num_subghosts_max_wave_speed_y[0]) +
                    (j - 1 + num_subghosts_max_wave_speed_y[1])*subghostcell_dims_max_wave_speed_y[0];
                
                const int idx_T_max_wave_speed_y = (i + num_subghosts_max_wave_speed_y[0]) +
                    (j + num_subghosts_max_wave_speed_y[1])*subghostcell_dims_max_wave_speed_y[0];
                
                const int idx_B_convective_flux_y = (i + num_subghosts_convective_flux_y[0]) +
                    (j - 1 + num_subghosts_convective_flux_y[1])*subghostcell_dims_convective_flux_y[0];
                
                const int idx_T_convective_flux_y = (i + num_subghosts_convective_flux_y[0]) +
                    (j + num_subghosts_convective_flux_y[1])*subghostcell_dims_convective_flux_y[0];
                
                const Real alpha_y = std::max(max_lambda_y[idx_B_max_wave_speed_y], max_lambda_y[idx_T_max_wave_speed_y]);
                
                // Compute the fluxes.
                for (int ei = 0; ei < d_num_eqn; ei++)
                {
                    // Compute the linear indices.
                    const int idx_B_conservative_var = (i + num_subghosts_conservative_var[ei][0]) +
                        (j - 1 + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0];
                    
                    const int idx_T_conservative_var = (i + num_subghosts_conservative_var[ei][0]) +
                        (j + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0];
                    
                    convective_flux->getPointer(1, ei)[idx_midpoint_y] = Real(1)/Real(2)*Real(dt)*(
                        F_y_node[ei][idx_B_convective_flux_y] + F_y_node[ei][idx_T_convective_flux_y] -
                            alpha_y*(Q[ei][idx_T_conservative_var] - Q[ei][idx_B_conservative_var]));
                }
                
                if (has_advection_eqn)
                {
                    // Compute the linear indices.
                    const int idx_B_density = (i + num_subghosts_density[0]) +
                        (j - 1 + num_subghosts_density[1])*subghostcell_dims_density[0];
                    
                    const int idx_T_density = (i + num_subghosts_density[0]) +
                        (j + num_subghosts_density[1])*subghostcell_dims_density[0];
                    
                    const int idx_B_pressure = (i + num_subghosts_pressure[0]) +
                        (j - 1 + num_subghosts_pressure[1])*subghostcell_dims_pressure[0];
                    
                    const int idx_T_pressure = (i + num_subghosts_pressure[0]) +
                        (j + num_subghosts_pressure[1])*subghostcell_dims_pressure[0];
                    
                    const int idx_B_velocity = (i + num_subghosts_velocity[0]) +
                        (j - 1 + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                    
                    const int idx_T_velocity = (i + num_subghosts_velocity[0]) +
                        (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0];
                    
                    velocity_intercell->getPointer(1, 0)[idx_midpoint_y] =
                        (rho[idx_B_density]*u[idx_B_velocity]*(-alpha_y - v[idx_B_velocity]) -
                            rho[idx_T_density]*u[idx_T_velocity]*(alpha_y - v[idx_T_velocity]))/
                                (rho[idx_B_density]*(-alpha_y - v[idx_B_velocity]) -
                                    rho[idx_T_density]*(alpha_y - v[idx_T_velocity]));
                    
                    velocity_intercell->getPointer(1, 1)[idx_midpoint_y] =
                        (p[idx_T_pressure] - p[idx_B_pressure] +
                            rho[idx_B_density]*v[idx_B_velocity]*(-alpha_y - v[idx_B_velocity]) -
                                rho[idx_T_density]*v[idx_T_velocity]*(alpha_y - v[idx_T_velocity]))/
                                    (rho[idx_B_density]*(-alpha_y - v[idx_B_velocity]) -
                                        rho[idx_T_density]*(alpha_y - v[idx_T_velocity]));
                }
            }
        }
        
        /*
         * Compute the source.
         */
        
        if (has_advection_eqn)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                if (eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = (i + num_subghosts_conservative_var[ei][0]) +
                                (j + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0];
                            
                            const int idx_cell_nghost = i + j*interior_dims[0];
                            
                            const int idx_midpoint_x_L = i +
                                j*(interior_dims[0] + 1);
                            
                            const int idx_midpoint_x_R = (i + 1) +
                                j*(interior_dims[0] + 1);
                            
                            const int idx_midpoint_y_B = i +
                                j*interior_dims[0];
                            
                            const int idx_midpoint_y_T = i +
                                (j + 1)*interior_dims[0];
                            
                            const Real& u_L = velocity_intercell->getPointer(0, 0)[idx_midpoint_x_L];
                            const Real& u_R = velocity_intercell->getPointer(0, 0)[idx_midpoint_x_R];
                            
                            const Real& v_B = velocity_intercell->getPointer(1, 1)[idx_midpoint_y_B];
                            const Real& v_T = velocity_intercell->getPointer(1, 1)[idx_midpoint_y_T];
                            
                            S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*((u_R - u_L)/Real(dx[0]) + (v_T - v_B)/Real(dx[1]));
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
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
         */
        
        d_flow_model->registerPatchWithDataContext(patch, data_context);
        
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        
        if (has_advection_eqn)
        {
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("DENSITY", d_num_conv_ghosts));
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRESSURE", d_num_conv_ghosts));
            num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("VELOCITY", d_num_conv_ghosts));
        }
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("MAX_WAVE_SPEED_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("MAX_WAVE_SPEED_Y", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("MAX_WAVE_SPEED_Z", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_X", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Y", d_num_conv_ghosts));
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("CONVECTIVE_FLUX_Z", d_num_conv_ghosts));
        
        d_flow_model->registerDerivedVariables(num_subghosts_of_data);
        
        d_flow_model->allocateMemoryForDerivedCellData();
        
        d_flow_model->computeDerivedCellData();
        
        /*
         * Get the pointer to the cell data inside the flow model.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > density;
        HAMERS_SHARED_PTR<pdat::CellData<Real> > pressure;
        HAMERS_SHARED_PTR<pdat::CellData<Real> > velocity;
        
        if (has_advection_eqn)
        {
            density = d_flow_model->getCellData("DENSITY");
            pressure = d_flow_model->getCellData("PRESSURE");
            velocity = d_flow_model->getCellData("VELOCITY");
        }
        
        HAMERS_SHARED_PTR<pdat::CellData<Real> > max_wave_speed_x = d_flow_model->getCellData("MAX_WAVE_SPEED_X");
        HAMERS_SHARED_PTR<pdat::CellData<Real> > max_wave_speed_y = d_flow_model->getCellData("MAX_WAVE_SPEED_Y");
        HAMERS_SHARED_PTR<pdat::CellData<Real> > max_wave_speed_z = d_flow_model->getCellData("MAX_WAVE_SPEED_Z");
        
        std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > convective_flux_node(3);
        convective_flux_node[0] = d_flow_model->getCellData("CONVECTIVE_FLUX_X");
        convective_flux_node[1] = d_flow_model->getCellData("CONVECTIVE_FLUX_Y");
        convective_flux_node[2] = d_flow_model->getCellData("CONVECTIVE_FLUX_Z");
        
        hier::IntVector num_subghosts_density(d_dim);
        hier::IntVector num_subghosts_pressure(d_dim);
        hier::IntVector num_subghosts_velocity(d_dim);
        
        hier::IntVector subghostcell_dims_density(d_dim);
        hier::IntVector subghostcell_dims_pressure(d_dim);
        hier::IntVector subghostcell_dims_velocity(d_dim);
        
        if (has_advection_eqn)
        {
            num_subghosts_density = density->getGhostCellWidth();
            subghostcell_dims_density = density->getGhostBox().numberCells();
            
            num_subghosts_pressure = pressure->getGhostCellWidth();
            subghostcell_dims_pressure = pressure->getGhostBox().numberCells();
            
            num_subghosts_velocity = velocity->getGhostCellWidth();
            subghostcell_dims_velocity = velocity->getGhostBox().numberCells();
        }
        
        hier::IntVector num_subghosts_max_wave_speed_x = max_wave_speed_x->getGhostCellWidth();
        hier::IntVector subghostcell_dims_max_wave_speed_x = max_wave_speed_x->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_max_wave_speed_y = max_wave_speed_y->getGhostCellWidth();
        hier::IntVector subghostcell_dims_max_wave_speed_y = max_wave_speed_y->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_max_wave_speed_z = max_wave_speed_z->getGhostCellWidth();
        hier::IntVector subghostcell_dims_max_wave_speed_z = max_wave_speed_z->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_x = convective_flux_node[0]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_x = convective_flux_node[0]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_y = convective_flux_node[1]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_y = convective_flux_node[1]->getGhostBox().numberCells();
        
        hier::IntVector num_subghosts_convective_flux_z = convective_flux_node[2]->getGhostCellWidth();
        hier::IntVector subghostcell_dims_convective_flux_z = convective_flux_node[2]->getGhostBox().numberCells();
        
        Real* rho = nullptr;
        Real* p = nullptr;
        Real* u = nullptr;
        Real* v = nullptr;
        Real* w = nullptr;
        
        if (has_advection_eqn)
        {
            rho = density->getPointer(0);
            p = pressure->getPointer(0);
            u = velocity->getPointer(0);
            v = velocity->getPointer(1);
            w = velocity->getPointer(2);
        }
        
        Real* max_lambda_x = max_wave_speed_x->getPointer(0);
        Real* max_lambda_y = max_wave_speed_y->getPointer(0);
        Real* max_lambda_z = max_wave_speed_z->getPointer(0);
        
        std::vector<Real*> F_x_node;
        std::vector<Real*> F_y_node;
        std::vector<Real*> F_z_node;
        F_x_node.reserve(d_num_eqn);
        F_y_node.reserve(d_num_eqn);
        F_z_node.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
            F_y_node.push_back(convective_flux_node[1]->getPointer(ei));
            F_z_node.push_back(convective_flux_node[2]->getPointer(ei));
        }
        
        /*
         * Get the pointers to the conservative variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
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
        
        for (int k = 0; k < interior_dims[2]; k++)
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0] + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_midpoint_x = i +
                        j*(interior_dims[0] + 1) +
                        k*(interior_dims[0] + 1)*interior_dims[1];
                    
                    const int idx_L_max_wave_speed_x = (i - 1 + num_subghosts_max_wave_speed_x[0]) +
                        (j + num_subghosts_max_wave_speed_x[1])*subghostcell_dims_max_wave_speed_x[0] +
                        (k + num_subghosts_max_wave_speed_x[2])*subghostcell_dims_max_wave_speed_x[0]*
                            subghostcell_dims_max_wave_speed_x[1];
                    
                    const int idx_R_max_wave_speed_x = (i + num_subghosts_max_wave_speed_x[0]) +
                        (j + num_subghosts_max_wave_speed_x[1])*subghostcell_dims_max_wave_speed_x[0] +
                        (k + num_subghosts_max_wave_speed_x[2])*subghostcell_dims_max_wave_speed_x[0]*
                            subghostcell_dims_max_wave_speed_x[1];
                    
                    const int idx_L_convective_flux_x = (i - 1 + num_subghosts_convective_flux_x[0]) +
                        (j + num_subghosts_convective_flux_x[1])*subghostcell_dims_convective_flux_x[0] +
                        (k + num_subghosts_convective_flux_x[2])*subghostcell_dims_convective_flux_x[0]*
                            subghostcell_dims_convective_flux_x[1];
                    
                    const int idx_R_convective_flux_x = (i + num_subghosts_convective_flux_x[0]) +
                        (j + num_subghosts_convective_flux_x[1])*subghostcell_dims_convective_flux_x[0] +
                        (k + num_subghosts_convective_flux_x[2])*subghostcell_dims_convective_flux_x[0]*
                            subghostcell_dims_convective_flux_x[1];
                    
                    const Real alpha_x = std::max(max_lambda_x[idx_L_max_wave_speed_x], max_lambda_x[idx_R_max_wave_speed_x]);
                    
                    // Compute the fluxes.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        // Compute the linear indices.
                        const int idx_L_conservative_var = (i - 1 + num_subghosts_conservative_var[ei][0]) +
                            (j + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0] +
                            (k + num_subghosts_conservative_var[ei][2])*subghostcell_dims_conservative_var[ei][0]*
                                subghostcell_dims_conservative_var[ei][1];
                        
                        const int idx_R_conservative_var = (i + num_subghosts_conservative_var[ei][0]) +
                            (j + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0] +
                            (k + num_subghosts_conservative_var[ei][2])*subghostcell_dims_conservative_var[ei][0]*
                                subghostcell_dims_conservative_var[ei][1];
                        
                        convective_flux->getPointer(0, ei)[idx_midpoint_x] = Real(1)/Real(2)*Real(dt)*(
                            F_x_node[ei][idx_L_convective_flux_x] + F_x_node[ei][idx_R_convective_flux_x] -
                                alpha_x*(Q[ei][idx_R_conservative_var] - Q[ei][idx_L_conservative_var]));
                    }
                    
                    if (has_advection_eqn)
                    {
                        // Compute the linear indices.
                        const int idx_L_density = (i - 1 + num_subghosts_density[0]) +
                            (j + num_subghosts_density[1])*subghostcell_dims_density[0] +
                            (k + num_subghosts_density[2])*subghostcell_dims_density[0]*
                                subghostcell_dims_density[1];
                        
                        const int idx_R_density = (i + num_subghosts_density[0]) +
                            (j + num_subghosts_density[1])*subghostcell_dims_density[0] +
                            (k + num_subghosts_density[2])*subghostcell_dims_density[0]*
                                subghostcell_dims_density[1];
                        
                        const int idx_L_pressure = (i - 1 + num_subghosts_pressure[0]) +
                            (j + num_subghosts_pressure[1])*subghostcell_dims_pressure[0] +
                            (k + num_subghosts_pressure[2])*subghostcell_dims_pressure[0]*
                                subghostcell_dims_pressure[1];
                        
                        const int idx_R_pressure = (i + num_subghosts_pressure[0]) +
                            (j + num_subghosts_pressure[1])*subghostcell_dims_pressure[0] +
                            (k + num_subghosts_pressure[2])*subghostcell_dims_pressure[0]*
                                subghostcell_dims_pressure[1];
                        
                        const int idx_L_velocity = (i - 1 + num_subghosts_velocity[0]) +
                            (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                            (k + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                subghostcell_dims_velocity[1];
                        
                        const int idx_R_velocity = (i + num_subghosts_velocity[0]) +
                            (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                            (k + num_subghosts_density[2])*subghostcell_dims_velocity[0]*
                                subghostcell_dims_density[1];
                        
                        velocity_intercell->getPointer(0, 0)[idx_midpoint_x] = (p[idx_R_pressure] - p[idx_L_pressure] +
                            rho[idx_L_density]*u[idx_L_velocity]*(-alpha_x - u[idx_L_velocity]) -
                                rho[idx_R_density]*u[idx_R_velocity]*(alpha_x - u[idx_R_velocity]))/
                                    (rho[idx_L_density]*(-alpha_x - u[idx_L_velocity]) -
                                        rho[idx_R_density]*(alpha_x - u[idx_R_velocity]));
                        
                        velocity_intercell->getPointer(0, 1)[idx_midpoint_x] =
                            (rho[idx_L_density]*v[idx_L_velocity]*(-alpha_x - u[idx_L_velocity]) -
                                rho[idx_R_density]*v[idx_R_velocity]*(alpha_x - u[idx_R_velocity]))/
                                    (rho[idx_L_density]*(-alpha_x - u[idx_L_velocity]) -
                                        rho[idx_R_density]*(alpha_x - u[idx_R_velocity]));
                        
                        velocity_intercell->getPointer(0, 2)[idx_midpoint_x] =
                            (rho[idx_L_density]*w[idx_L_velocity]*(-alpha_x - u[idx_L_velocity]) -
                                rho[idx_R_density]*w[idx_R_velocity]*(alpha_x - u[idx_R_velocity]))/
                                    (rho[idx_L_density]*(-alpha_x - u[idx_L_velocity]) -
                                        rho[idx_R_density]*(alpha_x - u[idx_R_velocity]));
                    }
                }
            }
        }
        
        /*
         * Compute the fluxes in the y direction.
         */
        
        for (int k = 0; k < interior_dims[2]; k++)
        {
            for (int j = 0; j < interior_dims[1] + 1; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute the linear indices.
                    const int idx_midpoint_y = i +
                        j*interior_dims[0] +
                        k*interior_dims[0]*(interior_dims[1] + 1);
                    
                    const int idx_B_max_wave_speed_y = (i + num_subghosts_max_wave_speed_y[0]) +
                        (j - 1 + num_subghosts_max_wave_speed_y[1])*subghostcell_dims_max_wave_speed_y[0] +
                        (k + num_subghosts_max_wave_speed_y[2])*subghostcell_dims_max_wave_speed_y[0]*
                            subghostcell_dims_max_wave_speed_y[1];
                    
                    const int idx_T_max_wave_speed_y = (i + num_subghosts_max_wave_speed_y[0]) +
                        (j + num_subghosts_max_wave_speed_y[1])*subghostcell_dims_max_wave_speed_y[0] +
                        (k + num_subghosts_max_wave_speed_y[2])*subghostcell_dims_max_wave_speed_y[0]*
                            subghostcell_dims_max_wave_speed_y[1];
                    
                    const int idx_B_convective_flux_y = (i + num_subghosts_convective_flux_y[0]) +
                        (j - 1 + num_subghosts_convective_flux_y[1])*subghostcell_dims_convective_flux_y[0] +
                        (k + num_subghosts_convective_flux_y[2])*subghostcell_dims_convective_flux_y[0]*
                            subghostcell_dims_convective_flux_y[1];
                    
                    const int idx_T_convective_flux_y = (i + num_subghosts_convective_flux_y[0]) +
                        (j + num_subghosts_convective_flux_y[1])*subghostcell_dims_convective_flux_y[0] +
                        (k + num_subghosts_convective_flux_y[2])*subghostcell_dims_convective_flux_y[0]*
                            subghostcell_dims_convective_flux_y[1];
                    
                    const Real alpha_y = std::max(max_lambda_y[idx_B_max_wave_speed_y], max_lambda_y[idx_T_max_wave_speed_y]);
                    
                    // Compute the fluxes.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        // Compute the linear indices.
                        const int idx_B_conservative_var = (i + num_subghosts_conservative_var[ei][0]) +
                        (j - 1 + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0] +
                        (k + num_subghosts_conservative_var[ei][2])*subghostcell_dims_conservative_var[ei][0]*
                            subghostcell_dims_conservative_var[ei][1];
                        
                        const int idx_T_conservative_var = (i + num_subghosts_conservative_var[ei][0]) +
                        (j + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0] +
                        (k + num_subghosts_conservative_var[ei][2])*subghostcell_dims_conservative_var[ei][0]*
                            subghostcell_dims_conservative_var[ei][1];
                        
                        convective_flux->getPointer(1, ei)[idx_midpoint_y] = Real(1)/Real(2)*Real(dt)*(
                            F_y_node[ei][idx_B_convective_flux_y] + F_y_node[ei][idx_T_convective_flux_y] -
                                alpha_y*(Q[ei][idx_T_conservative_var] - Q[ei][idx_B_conservative_var]));
                    }
                    
                    if (has_advection_eqn)
                    {
                        // Compute the linear indices.
                        const int idx_B_density = (i + num_subghosts_density[0]) +
                            (j - 1 + num_subghosts_density[1])*subghostcell_dims_density[0] +
                                (k + num_subghosts_density[2])*subghostcell_dims_density[0]*
                                    subghostcell_dims_density[1];
                        
                        const int idx_T_density = (i + num_subghosts_density[0]) +
                            (j + num_subghosts_density[1])*subghostcell_dims_density[0] +
                                (k + num_subghosts_density[2])*subghostcell_dims_density[0]*
                                    subghostcell_dims_density[1];
                        
                        const int idx_B_pressure = (i + num_subghosts_pressure[0]) +
                            (j - 1 + num_subghosts_pressure[1])*subghostcell_dims_pressure[0] +
                                (k + num_subghosts_pressure[2])*subghostcell_dims_pressure[0]*
                                    subghostcell_dims_pressure[1];
                        
                        const int idx_T_pressure = (i + num_subghosts_pressure[0]) +
                            (j + num_subghosts_pressure[1])*subghostcell_dims_pressure[0] +
                                (k + num_subghosts_pressure[2])*subghostcell_dims_pressure[0]*
                                    subghostcell_dims_pressure[1];
                        
                        const int idx_B_velocity = (i + num_subghosts_velocity[0]) +
                            (j - 1 + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                                (k + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                    subghostcell_dims_velocity[1];
                        
                        const int idx_T_velocity = (i + num_subghosts_velocity[0]) +
                            (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                                (k + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                    subghostcell_dims_velocity[1];
                        
                        velocity_intercell->getPointer(1, 0)[idx_midpoint_y] =
                            (rho[idx_B_density]*u[idx_B_velocity]*(-alpha_y - v[idx_B_velocity]) -
                                rho[idx_T_density]*u[idx_T_velocity]*(alpha_y - v[idx_T_velocity]))/
                                    (rho[idx_B_density]*(-alpha_y - v[idx_B_velocity]) -
                                        rho[idx_T_density]*(alpha_y - v[idx_T_velocity]));
                        
                        velocity_intercell->getPointer(1, 1)[idx_midpoint_y] =
                            (p[idx_T_pressure] - p[idx_B_pressure] +
                                rho[idx_B_density]*v[idx_B_velocity]*(-alpha_y - v[idx_B_velocity]) -
                                    rho[idx_T_density]*v[idx_T_velocity]*(alpha_y - v[idx_T_velocity]))/
                                        (rho[idx_B_density]*(-alpha_y - v[idx_B_velocity]) -
                                            rho[idx_T_density]*(alpha_y - v[idx_T_velocity]));
                        
                        velocity_intercell->getPointer(1, 2)[idx_midpoint_y] =
                            (rho[idx_B_density]*w[idx_B_velocity]*(-alpha_y - v[idx_B_velocity]) -
                                rho[idx_T_density]*w[idx_T_velocity]*(alpha_y - v[idx_T_velocity]))/
                                    (rho[idx_B_density]*(-alpha_y - v[idx_B_velocity]) -
                                        rho[idx_T_density]*(alpha_y - v[idx_T_velocity]));
                    }
                }
            }
        }
        
        /*
         * Compute the fluxes in the z direction.
         */
        
        for (int k = 0; k < interior_dims[2] + 1; k++)
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute the linear indices.
                    const int idx_midpoint_z = i +
                        j*interior_dims[0] +
                        k*interior_dims[0]*interior_dims[1];
                    
                    const int idx_B_max_wave_speed_z = (i + num_subghosts_max_wave_speed_z[0]) +
                        (j + num_subghosts_max_wave_speed_z[1])*subghostcell_dims_max_wave_speed_z[0] +
                        (k - 1 + num_subghosts_max_wave_speed_z[2])*subghostcell_dims_max_wave_speed_z[0]*
                            subghostcell_dims_max_wave_speed_z[1];
                        
                    const int idx_F_max_wave_speed_z = (i + num_subghosts_max_wave_speed_z[0]) +
                        (j + num_subghosts_max_wave_speed_z[1])*subghostcell_dims_max_wave_speed_z[0] +
                        (k + num_subghosts_max_wave_speed_z[2])*subghostcell_dims_max_wave_speed_z[0]*
                            subghostcell_dims_max_wave_speed_z[1];
                    
                    const int idx_B_convective_flux_z = (i + num_subghosts_convective_flux_z[0]) +
                        (j + num_subghosts_convective_flux_z[1])*subghostcell_dims_convective_flux_z[0] +
                        (k - 1 + num_subghosts_convective_flux_z[2])*subghostcell_dims_convective_flux_z[0]*
                            subghostcell_dims_convective_flux_z[1];
                        
                    const int idx_F_convective_flux_z = (i + num_subghosts_convective_flux_z[0]) +
                        (j + num_subghosts_convective_flux_z[1])*subghostcell_dims_convective_flux_z[0] +
                        (k + num_subghosts_convective_flux_z[2])*subghostcell_dims_convective_flux_z[0]*
                            subghostcell_dims_convective_flux_z[1];
                    
                    const Real alpha_z = std::max(max_lambda_z[idx_B_max_wave_speed_z], max_lambda_z[idx_F_max_wave_speed_z]);
                    
                    // Compute the fluxes.
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        const int idx_B_conservative_var = (i + num_subghosts_conservative_var[ei][0]) +
                            (j + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0] +
                            (k - 1 + num_subghosts_conservative_var[ei][2])*subghostcell_dims_conservative_var[ei][0]*
                                subghostcell_dims_conservative_var[ei][1];
                            
                        const int idx_F_conservative_var = (i + num_subghosts_conservative_var[ei][0]) +
                            (j + num_subghosts_conservative_var[ei][1])*subghostcell_dims_conservative_var[ei][0] +
                            (k + num_subghosts_conservative_var[ei][2])*subghostcell_dims_conservative_var[ei][0]*
                                subghostcell_dims_conservative_var[ei][1];
                        
                        convective_flux->getPointer(2, ei)[idx_midpoint_z] = Real(1)/Real(2)*Real(dt)*(
                            F_z_node[ei][idx_B_convective_flux_z] + F_z_node[ei][idx_F_convective_flux_z] -
                                alpha_z*(Q[ei][idx_F_conservative_var] - Q[ei][idx_B_conservative_var]));       
                    }
                    
                    if (has_advection_eqn)
                    {
                        // Compute the linear indices.
                        const int idx_B_density = (i + num_subghosts_density[0]) +
                            (j + num_subghosts_density[1])*subghostcell_dims_density[0] +
                            (k - 1 + num_subghosts_density[2])*subghostcell_dims_density[0]*
                                subghostcell_dims_density[1];
                            
                        const int idx_F_density = (i + num_subghosts_density[0]) +
                            (j + num_subghosts_density[1])*subghostcell_dims_density[0] +
                            (k + num_subghosts_density[2])*subghostcell_dims_density[0]*
                                subghostcell_dims_density[1];
                        
                        const int idx_B_pressure = (i + num_subghosts_pressure[0]) +
                            (j + num_subghosts_pressure[1])*subghostcell_dims_pressure[0] +
                            (k - 1 + num_subghosts_pressure[2])*subghostcell_dims_pressure[0]*
                                subghostcell_dims_pressure[1];
                            
                        const int idx_F_pressure = (i + num_subghosts_pressure[0]) +
                            (j + num_subghosts_pressure[1])*subghostcell_dims_pressure[0] +
                            (k + num_subghosts_pressure[2])*subghostcell_dims_pressure[0]*
                                subghostcell_dims_pressure[1];
                        
                        const int idx_B_velocity = (i + num_subghosts_velocity[0]) +
                            (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                            (k - 1 + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                subghostcell_dims_velocity[1];
                            
                        const int idx_F_velocity = (i + num_subghosts_velocity[0]) +
                            (j + num_subghosts_velocity[1])*subghostcell_dims_velocity[0] +
                                (k + num_subghosts_velocity[2])*subghostcell_dims_velocity[0]*
                                    subghostcell_dims_velocity[1];
                        
                        velocity_intercell->getPointer(2, 0)[idx_midpoint_z] =
                            (rho[idx_B_density]*u[idx_B_velocity]*(-alpha_z - w[idx_B_velocity]) -
                                rho[idx_F_density]*u[idx_F_velocity]*(alpha_z - w[idx_F_velocity]))/
                                    (rho[idx_B_density]*(-alpha_z - w[idx_B_velocity]) -
                                        rho[idx_F_density]*(alpha_z - w[idx_F_velocity]));
                        
                        velocity_intercell->getPointer(2, 1)[idx_midpoint_z] =
                            (rho[idx_B_density]*v[idx_B_velocity]*(-alpha_z - w[idx_B_velocity]) -
                                rho[idx_F_density]*v[idx_F_velocity]*(alpha_z - w[idx_F_velocity]))/
                                    (rho[idx_B_density]*(-alpha_z - w[idx_B_velocity]) -
                                        rho[idx_F_density]*(alpha_z - w[idx_F_density]));
                        
                        velocity_intercell->getPointer(2, 2)[idx_midpoint_z] =
                            (p[idx_F_pressure] - p[idx_B_pressure] +
                                rho[idx_B_density]*w[idx_B_velocity]*(-alpha_z - w[idx_B_velocity]) -
                                    rho[idx_F_density]*w[idx_F_velocity]*(alpha_z - w[idx_F_velocity]))/
                                        (rho[idx_B_density]*(-alpha_z - w[idx_B_velocity]) -
                                            rho[idx_F_density]*(alpha_z - w[idx_F_velocity]));
                    }
                }
            }
        }
        
        /*
         * Compute the source.
         */
        
        if (has_advection_eqn)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                if (eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    Real* S = source->getPointer(ei);
                    
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
                                
                                const int idx_midpoint_x_L = i +
                                    j*(interior_dims[0] + 1) +
                                    k*(interior_dims[0] + 1)*interior_dims[1];
                                
                                const int idx_midpoint_x_R = (i + 1) +
                                    j*(interior_dims[0] + 1) +
                                    k*(interior_dims[0] + 1)*interior_dims[1];
                                
                                const int idx_midpoint_y_B = i +
                                    j*interior_dims[0] +
                                    k*interior_dims[0]*(interior_dims[1] + 1);
                                
                                const int idx_midpoint_y_T = i +
                                    (j + 1)*interior_dims[0] +
                                    k*interior_dims[0]*(interior_dims[1] + 1);
                                
                                const int idx_midpoint_z_B = i +
                                    j*interior_dims[0] +
                                    k*interior_dims[0]*interior_dims[1];
                                
                                const int idx_midpoint_z_F = i +
                                    j*interior_dims[0] +
                                    (k + 1)*interior_dims[0]*interior_dims[1];
                                
                                const Real& u_L = velocity_intercell->getPointer(0, 0)[idx_midpoint_x_L];
                                const Real& u_R = velocity_intercell->getPointer(0, 0)[idx_midpoint_x_R];
                                
                                const Real& v_B = velocity_intercell->getPointer(1, 1)[idx_midpoint_y_B];
                                const Real& v_T = velocity_intercell->getPointer(1, 1)[idx_midpoint_y_T];
                                
                                const Real& w_B = velocity_intercell->getPointer(2, 2)[idx_midpoint_z_B];
                                const Real& w_F = velocity_intercell->getPointer(2, 2)[idx_midpoint_z_F];
                                
                                S[idx_cell_nghost] += Real(dt)*Q[ei][idx_cell_wghost]*(
                                    (u_R - u_L)/Real(dx[0]) + (v_T - v_B)/Real(dx[1]) + (w_F - w_B)/Real(dx[2]));
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
}
