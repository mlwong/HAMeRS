#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWENO-CU6-M2-LLF.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

ConvectiveFluxReconstructorWENO_CU6_M2_LLF::ConvectiveFluxReconstructorWENO_CU6_M2_LLF(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geom,
    const FLOW_MODEL& flow_model,
    const int& num_eqn,
    const int& num_species,
    const boost::shared_ptr<EquationOfState>& equation_of_state,
    const boost::shared_ptr<tbox::Database>& shock_capturing_scheme_db):
        ConvectiveFluxReconstructor(
            object_name,
            dim,
            grid_geom,
            flow_model,
            num_eqn,
            num_species,
            equation_of_state,
            shock_capturing_scheme_db)
{
    d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*3;
    
    d_constant_C       = d_shock_capturing_scheme_db->getDoubleWithDefault("constant_C", 1000.0);
    d_constant_C       = d_shock_capturing_scheme_db->getDoubleWithDefault("d_constant_C", d_constant_C);
    
    d_constant_q       = d_shock_capturing_scheme_db->getIntegerWithDefault("constant_q", 4);
    d_constant_q       = d_shock_capturing_scheme_db->getIntegerWithDefault("d_constant_q", d_constant_q);
    
    d_constant_epsilon = d_shock_capturing_scheme_db->getDoubleWithDefault("constant_epsilon", 1.0e-8);
    d_constant_epsilon = d_shock_capturing_scheme_db->getDoubleWithDefault("d_constant_epsilon", d_constant_epsilon);
    
    d_constant_Chi     = d_shock_capturing_scheme_db->getDoubleWithDefault("constant_Chi", 1.0e8);
    d_constant_Chi     = d_shock_capturing_scheme_db->getDoubleWithDefault("d_constant_Chi", d_constant_Chi);
    
    d_weights_d.push_back(1.0/20);
    d_weights_d.push_back(9.0/20);
    d_weights_d.push_back(9.0/20);
    d_weights_d.push_back(1.0/20);
    
    d_weights_c.resize(boost::extents[4][3]);
    d_weights_c[0][0] = 1.0/3;
    d_weights_c[0][1] = -7.0/6;
    d_weights_c[0][2] = 11.0/6;
    d_weights_c[1][0] = -1.0/6;
    d_weights_c[1][1] = 5.0/6;
    d_weights_c[1][2] = 1.0/3;
    d_weights_c[2][0] = 1.0/3;
    d_weights_c[2][1] = 5.0/6;
    d_weights_c[2][2] = -1.0/6;
    d_weights_c[3][0] = 11.0/6;
    d_weights_c[3][1] = -7.0/6;
    d_weights_c[3][2] = 1.0/3;
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorWENO_CU6_M2_LLF::printClassData(
    std::ostream& os) const
{
    os << "\nConvectiveFluxReconstructorWENO_CU6_M2_LLF::printClassData..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorWENO_CU6_M2_LLF: this = "
       << (ConvectiveFluxReconstructorWENO_CU6_M2_LLF *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_set_variables = "
       << d_set_variables
       << std::endl;
    os << "d_constant_C = "
       << d_constant_C
       << std::endl;
    os << "d_constant_q = "
       << d_constant_q
       << std::endl;
    
    os << std::endl;
    
    os << "End of ConvectiveFluxReconstructorWENO_CU6_M2_LLF::printClassData"
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorWENO_CU6_M2_LLF::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putString("d_shock_capturing_scheme", "WENO_CU6_M2_LLF");
    
    restart_db->putDouble("d_constant_C", d_constant_C);
    restart_db->putInteger("d_constant_q", d_constant_q);
    restart_db->putDouble("d_constant_epsilon", d_constant_epsilon);
    restart_db->putDouble("d_constant_Chi", d_constant_Chi);
}


/*
 * Compute the convective fluxes and sources due to hyperbolization
 * of the equations.
 */
void
ConvectiveFluxReconstructorWENO_CU6_M2_LLF::computeConvectiveFluxesAndSources(
    hier::Patch& patch,
    const double time,
    const double dt,
    const int RK_step_number,
    const boost::shared_ptr<hier::VariableContext> data_context)
{
    NULL_USE(RK_step_number);
    
    if (d_set_variables == true)
    {
        // Get the dimensions of box that covers the interior of patch.
        hier::Box dummy_box = patch.getBox();
        const hier::Box interior_box = dummy_box;
        const hier::IntVector interior_dims = interior_box.numberCells();
        
        // Get the dimensions of box that covers interior of patch plus
        // ghost cells.
        dummy_box.grow(d_num_ghosts);
        const hier::Box ghost_box = dummy_box;
        const hier::IntVector ghostcell_dims = ghost_box.numberCells();
        
        switch (d_flow_model)
        {
            case SINGLE_SPECIES:
            {
                if (d_equation_of_state->getLabel() == IDEAL_GAS)
                {
                    // Get the cell data of the time-dependent variables.
                    boost::shared_ptr<pdat::CellData<double> > density(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_density, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > momentum(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_momentum, data_context)));
                    
                    boost::shared_ptr<pdat::CellData<double> > total_energy(
                        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                            patch.getPatchData(d_total_energy, data_context)));
                    
                    boost::shared_ptr<pdat::FaceData<double> > convective_flux(
                        BOOST_CAST<pdat::FaceData<double>, hier::PatchData>(
                            patch.getPatchData(d_convective_flux, data_context)));
                    
#ifdef DEBUG_CHECK_ASSERTIONS
                    TBOX_ASSERT(density);
                    TBOX_ASSERT(momentum);
                    TBOX_ASSERT(total_energy);
                    TBOX_ASSERT(convective_flux);
                    
                    TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
                    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
                    
                    // Allocate temporary patch data.
                    boost::shared_ptr<pdat::CellData<double> > velocity(
                        new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > pressure(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > sound_speed(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > total_specific_enthalpy(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                    
                    std::vector<boost::shared_ptr<pdat::CellData<double> > > convective_flux_node;
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        convective_flux_node.push_back(boost::make_shared<pdat::CellData<double> >(
                            interior_box, d_num_eqn, d_num_ghosts));
                    }
                    
                    std::vector<boost::shared_ptr<pdat::CellData<double> > > wave_speed;
                    for (int di = 0; di < d_dim.getValue(); di++)
                    {
                        wave_speed.push_back(boost::make_shared<pdat::CellData<double> >(
                            interior_box, d_num_eqn, d_num_ghosts));
                    }
                    
                    boost::shared_ptr<pdat::FaceData<double> > projection_matrix(
                    new pdat::FaceData<double>(
                        interior_box,
                        d_num_eqn*d_num_eqn,
                        hier::IntVector::getZero(d_dim)));
                    
                    boost::shared_ptr<pdat::FaceData<double> > projection_matrix_inv(
                        new pdat::FaceData<double>(
                            interior_box,
                            d_num_eqn*d_num_eqn,
                            hier::IntVector::getZero(d_dim)));
                    
                    if (d_dim == tbox::Dimension(1))
                    {
                        // Get the arrays of time-dependent variables.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* E     = total_energy->getPointer(0);
                        
                        // Get the arrays of temporary patch data.
                        double* u     = velocity->getPointer(0);
                        double* p     = pressure->getPointer(0);
                        double* c     = sound_speed->getPointer(0);
                        double* H     = total_specific_enthalpy->getPointer(0);
                        std::vector<double*> F_x_node;
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
                        }
                        std::vector<double*> lambda_x;
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            lambda_x.push_back(wave_speed[0]->getPointer(ei));
                        }
                        
                        // Compute the field of velocities, pressure, sound speed, total specific enthalpy, wave speeds
                        // and fluxes.
                        for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                        {
                            // Compute index into linear data array.
                            const int idx = i + d_num_ghosts[0];
                            
                            u[idx] = rho_u[idx]/rho[idx];
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&rho_u[idx]);
                            
                            p[idx] = d_equation_of_state->
                                getPressure(
                                    &rho[idx],
                                    m_ptr,
                                    &E[idx]);
                            
                            c[idx] = d_equation_of_state->
                                getSoundSpeedWithPressure(
                                    &rho[idx],
                                    &p[idx]);
                            
                            H[idx] = (E[idx] + p[idx])/rho[idx];
                            
                            F_x_node[0][idx] = rho_u[idx];
                            F_x_node[1][idx] = rho_u[idx]*u[idx] + p[idx];
                            F_x_node[2][idx] = u[idx]*(E[idx] + p[idx]);
                            
                            lambda_x[0][idx] = u[idx] - c[idx];
                            lambda_x[1][idx] = u[idx];
                            lambda_x[2][idx] = u[idx] + c[idx];
                        }
                        
                        /*
                         * Compute the projection matrix and its inverse at the face normal to the
                         * x direction.
                         */
                        for (int i = 0; i < interior_dims[0] + 1; i++)
                        {
                            // Compute the indices of left cell, right cell and face of
                            // projection matrix.
                            const int idx_cell_L = i - 1 + d_num_ghosts[0];
                            
                            const int idx_cell_R = i + d_num_ghosts[0];
                            
                            const int idx_face_x = i;
                            
                            // Get left and right quantities.
                            const double& u_L = u[idx_cell_L];
                            const double& u_R = u[idx_cell_R];
                            
                            const double& c_L = c[idx_cell_L];
                            const double& c_R = c[idx_cell_R];
                            
                            const double& H_L = H[idx_cell_L];
                            const double& H_R = H[idx_cell_R];
                            
                            // Compute simply-averaged quantities.
                            const double u_average = 0.5*(u_L + u_R);
                            const double c_average = 0.5*(c_L + c_R);
                            const double H_average = 0.5*(H_L + H_R);
                            
                            boost::multi_array<double*, 2> R_x_intercell(
                                boost::extents[d_num_eqn][d_num_eqn]);
                            
                            boost::multi_array<double*, 2> R_x_inv_intercell(
                                boost::extents[d_num_eqn][d_num_eqn]);
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                for (int ej = 0; ej < d_num_eqn; ej++)
                                {
                                    R_x_intercell[ei][ej] =
                                        &(projection_matrix_inv->getPointer(0, ei + ej*d_num_eqn)[idx_face_x]);
                                    
                                    R_x_inv_intercell[ei][ej] =
                                        &(projection_matrix->getPointer(0, ei + ej*d_num_eqn)[idx_face_x]);
                                }
                            }
                            
                            *R_x_intercell[0][0] = 1.0;
                            *R_x_intercell[0][1] = 1.0;
                            *R_x_intercell[0][2] = 1.0;
                            *R_x_intercell[1][0] = u_average - c_average;
                            *R_x_intercell[1][1] = u_average;
                            *R_x_intercell[1][2] = u_average + c_average;
                            *R_x_intercell[2][0] = H_average - u_average*c_average;
                            *R_x_intercell[2][1] = 0.5*u_average*u_average;
                            *R_x_intercell[2][2] = H_average + u_average*c_average;
                            
                            const double gamma = d_equation_of_state->
                                getSpeciesThermodynamicProperty(
                                    "gamma",
                                    0);
                            const double b_1 = (gamma - 1.0)/(c_average*c_average);
                            const double b_2 = 0.5*u_average*u_average*b_1;
                            
                            *R_x_inv_intercell[0][0] = 0.5*(b_2 + u_average/c_average);
                            *R_x_inv_intercell[0][1] = -0.5*(b_1*u_average + 1.0/c_average);
                            *R_x_inv_intercell[0][2] = 0.5*b_1;
                            *R_x_inv_intercell[1][0] = 1.0 - b_2;
                            *R_x_inv_intercell[1][1] = b_1*u_average;
                            *R_x_inv_intercell[1][2] = -b_1;
                            *R_x_inv_intercell[2][0] = 0.5*(b_2 - u_average/c_average);
                            *R_x_inv_intercell[2][1] = 0.5*(-b_1*u_average + 1.0/c_average);
                            *R_x_inv_intercell[2][2] = 0.5*b_1;
                        }
                        
                        // Compute the fluxes in the x direction.
                        for (int i = 0; i < interior_dims[0] + 1; i++)
                        {
                            // Compute the indices of left cell, right cell and face.
                            const int idx_cell_L = i - 1 + d_num_ghosts[0];
                            
                            const int idx_cell_R = i + d_num_ghosts[0];
                        
                            const int idx_face_x = i;
                            
                            /*
                             * Compute the maximum wave speed of each equation locally.
                             */
                            
                            std::vector<double> alpha_x;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                alpha_x.push_back(fmax(fabs(lambda_x[ei][idx_cell_L]),
                                    fabs(lambda_x[ei][idx_cell_R])));
                            }
                            
                            boost::multi_array<double, 2> G_x_plus_array(
                                boost::extents[6][d_num_eqn]);
                            
                            boost::multi_array<double, 2> G_x_minus_array(
                                boost::extents[6][d_num_eqn]);
                            
                            /*
                             * Project fluxes onto characteristic fields.
                             */
                            
                            boost::multi_array<const double*, 2> R_x_inv_intercell(
                                boost::extents[d_num_eqn][d_num_eqn]);
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                for (int ej = 0; ej < d_num_eqn; ej++)
                                {
                                    R_x_inv_intercell[ei][ej] =
                                        &(projection_matrix->getPointer(0, ei + ej*d_num_eqn)[idx_face_x]);
                                }
                            }
                            
                            for (int m = 0; m < 6; m++)
                            {
                                const int idx_cell = i - 3 + m + d_num_ghosts[0];
                                
                                std::vector<const double*> F_x_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_x_ptr.push_back(&F_x_node[ei][idx_cell]);
                                }
                                
                                std::vector<double> G_x(d_num_eqn, 0.0);
                                std::vector<double*> G_x_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    G_x_ptr.push_back(&G_x[ei]);
                                }
                                
                                projectVariablesToCharacteristicFields(
                                    G_x_ptr,
                                    F_x_ptr,
                                    R_x_inv_intercell);
                                
                                std::vector<const double*> Q_ptr;
                                Q_ptr.push_back(&rho[idx_cell]);
                                Q_ptr.push_back(&rho_u[idx_cell]);
                                Q_ptr.push_back(&E[idx_cell]);
                                
                                std::vector<double> V(d_num_eqn, 0.0);
                                std::vector<double*> V_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    V_ptr.push_back(&V[ei]);
                                }
                                
                                projectVariablesToCharacteristicFields(
                                    V_ptr,
                                    Q_ptr,
                                    R_x_inv_intercell);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    G_x_plus_array[m][ei] = 0.5*(G_x[ei] +
                                        alpha_x[ei]*V[ei]);
                                    
                                    G_x_minus_array[m][ei] = 0.5*(G_x[ei] -
                                        alpha_x[ei]*V[ei]);
                                }
                            }
                            
                            /*
                             * Do WENO reconstruction on characteristic variables to get G_x_plus
                             * and G_x_minus
                             */
                            
                            std::vector<double> G_x_plus;
                            std::vector<double> G_x_minus;
                            
                            performWENOReconstruction(
                                G_x_plus,
                                G_x_minus,
                                G_x_plus_array,
                                G_x_minus_array,
                                X_DIRECTION);
                            
                            /*
                             * Project characteristic variables back to physical fields.
                             */
                            
                            boost::multi_array<const double*, 2> R_x_intercell(
                                boost::extents[d_num_eqn][d_num_eqn]);
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                for (int ej = 0; ej < d_num_eqn; ej++)
                                {
                                    R_x_intercell[ei][ej] =
                                        &(projection_matrix_inv->getPointer(0, ei + ej*d_num_eqn)[idx_face_x]);
                                }
                            }
                            
                            std::vector<double> F_x_plus(d_num_eqn, 0.0);
                            std::vector<double> F_x_minus(d_num_eqn, 0.0);
                            
                            std::vector<double*> F_x_plus_ptr;
                            std::vector<double*> F_x_minus_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_x_plus_ptr.push_back(&F_x_plus[ei]);
                                F_x_minus_ptr.push_back(&F_x_minus[ei]);
                            }
                            
                            std::vector<const double*> G_x_plus_ptr;
                            std::vector<const double*> G_x_minus_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                G_x_plus_ptr.push_back(&G_x_plus[ei]);
                                G_x_minus_ptr.push_back(&G_x_minus[ei]);
                            }
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                F_x_plus_ptr,
                                G_x_plus_ptr,
                                R_x_intercell);
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                F_x_minus_ptr,
                                G_x_minus_ptr,
                                R_x_intercell);
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                double* F_x = convective_flux->getPointer(0, ei);
                                F_x[idx_face_x] = dt*(F_x_plus[ei] + F_x_minus[ei]);
                            }
                        }
                    }
                    else if (d_dim == tbox::Dimension(2))
                    {
                        // Get the arrays of time-dependent variables.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* E     = total_energy->getPointer(0);
                        
                        // Get the arrays of temporary patch data.
                        double* u     = velocity->getPointer(0);
                        double* v     = velocity->getPointer(1);
                        double* p     = pressure->getPointer(0);
                        double* c     = sound_speed->getPointer(0);
                        double* H     = total_specific_enthalpy->getPointer(0);
                        std::vector<double*> F_x_node;
                        std::vector<double*> F_y_node;
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
                            F_y_node.push_back(convective_flux_node[1]->getPointer(ei));
                        }
                        std::vector<double*> lambda_x;
                        std::vector<double*> lambda_y;
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            lambda_x.push_back(wave_speed[0]->getPointer(ei));
                            lambda_y.push_back(wave_speed[1]->getPointer(ei));
                        }
                        
                        // Compute the field of velocities, pressure, sound speed, total specific enthalpy, wave speeds
                        // and fluxes.
                        for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                        {
                            for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                            {
                                // Compute index into linear data array.
                                const int idx = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                u[idx] = rho_u[idx]/rho[idx];
                                v[idx] = rho_v[idx]/rho[idx];
                                
                                std::vector<const double*> m_ptr;
                                m_ptr.push_back(&rho_u[idx]);
                                m_ptr.push_back(&rho_v[idx]);
                                
                                p[idx] = d_equation_of_state->
                                    getPressure(
                                        &rho[idx],
                                        m_ptr,
                                        &E[idx]);
                                
                                c[idx] = d_equation_of_state->
                                    getSoundSpeedWithPressure(
                                        &rho[idx],
                                        &p[idx]);
                                
                                H[idx] = (E[idx] + p[idx])/rho[idx];
                                
                                F_x_node[0][idx] = rho_u[idx];
                                F_x_node[1][idx] = rho_u[idx]*u[idx] + p[idx];
                                F_x_node[2][idx] = rho_u[idx]*v[idx];
                                F_x_node[3][idx] = u[idx]*(E[idx] + p[idx]);
                                
                                F_y_node[0][idx] = rho_v[idx];
                                F_y_node[1][idx] = rho_v[idx]*u[idx];
                                F_y_node[2][idx] = rho_v[idx]*v[idx] + p[idx];
                                F_y_node[3][idx] = v[idx]*(E[idx] + p[idx]);
                                
                                lambda_x[0][idx] = u[idx] - c[idx];
                                lambda_x[1][idx] = u[idx];
                                lambda_x[2][idx] = u[idx];
                                lambda_x[3][idx] = u[idx] + c[idx];
                                
                                lambda_y[0][idx] = v[idx] - c[idx];
                                lambda_y[1][idx] = v[idx];
                                lambda_y[2][idx] = v[idx];
                                lambda_y[3][idx] = v[idx] + c[idx];
                            }
                        }
                        
                        /*
                         * Compute the projection matrix and its inverse at the face normal to the
                         * x direction.
                         */
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0] + 1; i++)
                            {
                                // Compute the indices of left cell, right cell and face of
                                // projection matrix.
                                const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_R = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_face_x = i +
                                    j*(interior_dims[0] + 1);
                                
                                // Get left and right quantities.
                                const double& u_L = u[idx_cell_L];
                                const double& u_R = u[idx_cell_R];
                                
                                const double& v_L = v[idx_cell_L];
                                const double& v_R = v[idx_cell_R];
                                
                                const double& c_L = c[idx_cell_L];
                                const double& c_R = c[idx_cell_R];
                                
                                const double& H_L = H[idx_cell_L];
                                const double& H_R = H[idx_cell_R];
                                
                                // Compute simply-averaged quantities.
                                const double u_average = 0.5*(u_L + u_R);
                                const double v_average = 0.5*(v_L + v_R);
                                const double c_average = 0.5*(c_L + c_R);
                                const double H_average = 0.5*(H_L + H_R);
                                
                                boost::multi_array<double*, 2> R_x_intercell(
                                    boost::extents[d_num_eqn][d_num_eqn]);
                                
                                boost::multi_array<double*, 2> R_x_inv_intercell(
                                    boost::extents[d_num_eqn][d_num_eqn]);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    for (int ej = 0; ej < d_num_eqn; ej++)
                                    {
                                        R_x_intercell[ei][ej] =
                                            &(projection_matrix_inv->getPointer(0, ei + ej*d_num_eqn)[idx_face_x]);
                                        
                                        R_x_inv_intercell[ei][ej] =
                                            &(projection_matrix->getPointer(0, ei + ej*d_num_eqn)[idx_face_x]);
                                    }
                                }
                                
                                *R_x_intercell[0][0] = 1.0;
                                *R_x_intercell[0][1] = 0.0;
                                *R_x_intercell[0][2] = 1.0;
                                *R_x_intercell[0][3] = 1.0;
                                *R_x_intercell[1][0] = u_average - c_average;
                                *R_x_intercell[1][1] = 0.0;
                                *R_x_intercell[1][2] = u_average;
                                *R_x_intercell[1][3] = u_average + c_average;
                                *R_x_intercell[2][0] = v_average;
                                *R_x_intercell[2][1] = 1.0;
                                *R_x_intercell[2][2] = v_average;
                                *R_x_intercell[2][3] = v_average;
                                *R_x_intercell[3][0] = H_average - u_average*c_average;
                                *R_x_intercell[3][1] = v_average;
                                *R_x_intercell[3][2] = 0.5*(u_average*u_average + v_average*v_average);
                                *R_x_intercell[3][3] = H_average + u_average*c_average;
                                
                                const double gamma = d_equation_of_state->
                                    getSpeciesThermodynamicProperty(
                                        "gamma",
                                        0);
                                const double b_1 = (gamma - 1.0)/(c_average*c_average);
                                const double b_2 = 0.5*(u_average*u_average + v_average*v_average)*b_1;
                                
                                *R_x_inv_intercell[0][0] = 0.5*(b_2 + u_average/c_average);
                                *R_x_inv_intercell[0][1] = -0.5*(b_1*u_average + 1.0/c_average);
                                *R_x_inv_intercell[0][2] = -0.5*b_1*v_average;
                                *R_x_inv_intercell[0][3] = 0.5*b_1;
                                *R_x_inv_intercell[1][0] = -v_average;
                                *R_x_inv_intercell[1][1] = 0.0;
                                *R_x_inv_intercell[1][2] = 1.0;
                                *R_x_inv_intercell[1][3] = 0.0;
                                *R_x_inv_intercell[2][0] = 1.0 - b_2;
                                *R_x_inv_intercell[2][1] = b_1*u_average;
                                *R_x_inv_intercell[2][2] = b_1*v_average;
                                *R_x_inv_intercell[2][3] = -b_1;
                                *R_x_inv_intercell[3][0] = 0.5*(b_2 - u_average/c_average);
                                *R_x_inv_intercell[3][1] = 0.5*(-b_1*u_average + 1.0/c_average);
                                *R_x_inv_intercell[3][2] = -0.5*b_1*v_average;
                                *R_x_inv_intercell[3][3] = 0.5*b_1;
                            }
                        }
                        
                        /*
                         * Compute the projection matrix and its inverse at the face normal to the
                         * y direction.
                         */
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int j = 0; j < interior_dims[1] + 1; j++)
                            {
                                // Compute the indices of bottom cell, top cell and face of
                                // projection matrix.
                                const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_T = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_face_y = j +
                                    i*(interior_dims[1] + 1);
                                
                                // Get bottom and top quantities.
                                const double& u_B = u[idx_cell_B];
                                const double& u_T = u[idx_cell_T];
                                
                                const double& v_B = v[idx_cell_B];
                                const double& v_T = v[idx_cell_T];
                                
                                const double& c_B = c[idx_cell_B];
                                const double& c_T = c[idx_cell_T];
                                
                                const double& H_B = H[idx_cell_B];
                                const double& H_T = H[idx_cell_T];
                                
                                // Compute simply-averaged quantities.
                                const double u_average = 0.5*(u_B + u_T);
                                const double v_average = 0.5*(v_B + v_T);
                                const double c_average = 0.5*(c_B + c_T);
                                const double H_average = 0.5*(H_B + H_T);
                                
                                boost::multi_array<double*, 2> R_y_intercell(
                                    boost::extents[d_num_eqn][d_num_eqn]);
                                
                                boost::multi_array<double*, 2> R_y_inv_intercell(
                                    boost::extents[d_num_eqn][d_num_eqn]);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    for (int ej = 0; ej < d_num_eqn; ej++)
                                    {
                                        R_y_intercell[ei][ej] =
                                            &(projection_matrix_inv->getPointer(1, ei + ej*d_num_eqn)[idx_face_y]);
                                        
                                        R_y_inv_intercell[ei][ej] =
                                            &(projection_matrix->getPointer(1, ei + ej*d_num_eqn)[idx_face_y]);
                                    }
                                }
                                
                                *R_y_intercell[0][0] = 1.0;
                                *R_y_intercell[0][1] = 0.0;
                                *R_y_intercell[0][2] = 1.0;
                                *R_y_intercell[0][3] = 1.0;
                                *R_y_intercell[1][0] = u_average;
                                *R_y_intercell[1][1] = 1.0;
                                *R_y_intercell[1][2] = u_average;
                                *R_y_intercell[1][3] = u_average;
                                *R_y_intercell[2][0] = v_average - c_average;
                                *R_y_intercell[2][1] = 0.0;
                                *R_y_intercell[2][2] = v_average;
                                *R_y_intercell[2][3] = v_average + c_average;
                                *R_y_intercell[3][0] = H_average - v_average*c_average;
                                *R_y_intercell[3][1] = u_average;
                                *R_y_intercell[3][2] = 0.5*(u_average*u_average + v_average*v_average);
                                *R_y_intercell[3][3] = H_average + v_average*c_average;
                                
                                const double gamma = d_equation_of_state->
                                    getSpeciesThermodynamicProperty(
                                        "gamma",
                                        0);
                                const double b_1 = (gamma - 1.0)/(c_average*c_average);
                                const double b_2 = 0.5*(u_average*u_average + v_average*v_average)*b_1;
                                
                                *R_y_inv_intercell[0][0] = 0.5*(b_2 + v_average/c_average);
                                *R_y_inv_intercell[0][1] = -0.5*b_1*u_average;
                                *R_y_inv_intercell[0][2] = -0.5*(b_1*v_average + 1.0/c_average);
                                *R_y_inv_intercell[0][3] = 0.5*b_1;
                                *R_y_inv_intercell[1][0] = -u_average;
                                *R_y_inv_intercell[1][1] = 1.0;
                                *R_y_inv_intercell[1][2] = 0.0;
                                *R_y_inv_intercell[1][3] = 0.0;
                                *R_y_inv_intercell[2][0] = 1.0 - b_2;
                                *R_y_inv_intercell[2][1] = b_1*u_average;
                                *R_y_inv_intercell[2][2] = b_1*v_average;
                                *R_y_inv_intercell[2][3] = -b_1;
                                *R_y_inv_intercell[3][0] = 0.5*(b_2 - v_average/c_average);
                                *R_y_inv_intercell[3][1] = -0.5*b_1*u_average;
                                *R_y_inv_intercell[3][2] = 0.5*(-b_1*v_average + 1.0/c_average);
                                *R_y_inv_intercell[3][3] = 0.5*b_1;
                            }
                        }
                        
                        // Compute the fluxes in the x direction.
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0] + 1; i++)
                            {
                                // Compute the indices of left cell, right cell and face.
                                const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_R = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                                const int idx_face_x = i +
                                    j*(interior_dims[0] + 1);
                                
                                /*
                                 * Compute the maximum wave speed of each equation locally.
                                 */
                                
                                std::vector<double> alpha_x;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    alpha_x.push_back(fmax(fabs(lambda_x[ei][idx_cell_L]),
                                        fabs(lambda_x[ei][idx_cell_R])));
                                }
                                
                                boost::multi_array<double, 2> G_x_plus_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                boost::multi_array<double, 2> G_x_minus_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                /*
                                 * Project fluxes onto characteristic fields.
                                 */
                                
                                boost::multi_array<const double*, 2> R_x_inv_intercell(
                                    boost::extents[d_num_eqn][d_num_eqn]);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    for (int ej = 0; ej < d_num_eqn; ej++)
                                    {
                                        R_x_inv_intercell[ei][ej] =
                                            &(projection_matrix->getPointer(0, ei + ej*d_num_eqn)[idx_face_x]);
                                    }
                                }
                                
                                for (int m = 0; m < 6; m++)
                                {
                                    const int idx_cell = (i - 3 + m + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0];
                                    
                                    std::vector<const double*> F_x_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        F_x_ptr.push_back(&F_x_node[ei][idx_cell]);
                                    }
                                    
                                    std::vector<double> G_x(d_num_eqn, 0.0);
                                    std::vector<double*> G_x_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        G_x_ptr.push_back(&G_x[ei]);
                                    }
                                    
                                    projectVariablesToCharacteristicFields(
                                        G_x_ptr,
                                        F_x_ptr,
                                        R_x_inv_intercell);
                                    
                                    std::vector<const double*> Q_ptr;
                                    Q_ptr.push_back(&rho[idx_cell]);
                                    Q_ptr.push_back(&rho_u[idx_cell]);
                                    Q_ptr.push_back(&rho_v[idx_cell]);
                                    Q_ptr.push_back(&E[idx_cell]);
                                    
                                    std::vector<double> V(d_num_eqn, 0.0);
                                    std::vector<double*> V_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        V_ptr.push_back(&V[ei]);
                                    }
                                    
                                    projectVariablesToCharacteristicFields(
                                        V_ptr,
                                        Q_ptr,
                                        R_x_inv_intercell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        G_x_plus_array[m][ei] = 0.5*(G_x[ei] +
                                            alpha_x[ei]*V[ei]);
                                        
                                        G_x_minus_array[m][ei] = 0.5*(G_x[ei] -
                                            alpha_x[ei]*V[ei]);
                                    }
                                }
                                
                                /*
                                 * Do WENO reconstruction on characteristic variables to get G_x_plus
                                 * and G_x_minus
                                 */
                                
                                std::vector<double> G_x_plus;
                                std::vector<double> G_x_minus;
                                
                                performWENOReconstruction(
                                    G_x_plus,
                                    G_x_minus,
                                    G_x_plus_array,
                                    G_x_minus_array,
                                    X_DIRECTION);
                                
                                /*
                                 * Project characteristic variables back to physical fields.
                                 */
                                
                                boost::multi_array<const double*, 2> R_x_intercell(
                                    boost::extents[d_num_eqn][d_num_eqn]);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    for (int ej = 0; ej < d_num_eqn; ej++)
                                    {
                                        R_x_intercell[ei][ej] =
                                            &(projection_matrix_inv->getPointer(0, ei + ej*d_num_eqn)[idx_face_x]);
                                    }
                                }
                                
                                std::vector<double> F_x_plus(d_num_eqn, 0.0);
                                std::vector<double> F_x_minus(d_num_eqn, 0.0);
                                
                                std::vector<double*> F_x_plus_ptr;
                                std::vector<double*> F_x_minus_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_x_plus_ptr.push_back(&F_x_plus[ei]);
                                    F_x_minus_ptr.push_back(&F_x_minus[ei]);
                                }
                                
                                std::vector<const double*> G_x_plus_ptr;
                                std::vector<const double*> G_x_minus_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    G_x_plus_ptr.push_back(&G_x_plus[ei]);
                                    G_x_minus_ptr.push_back(&G_x_minus[ei]);
                                }
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    F_x_plus_ptr,
                                    G_x_plus_ptr,
                                    R_x_intercell);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    F_x_minus_ptr,
                                    G_x_minus_ptr,
                                    R_x_intercell);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double* F_x = convective_flux->getPointer(0, ei);
                                    F_x[idx_face_x] = dt*(F_x_plus[ei] + F_x_minus[ei]);
                                }
                            }
                        }
                        
                        // Compute the fluxes in the y direction.
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int j = 0; j < interior_dims[1] + 1; j++)
                            {
                                // Compute the indices of bottom cell, top cell and face.
                                const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_T = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_face_y = j +
                                    i*(interior_dims[1] + 1);
                                
                                /*
                                 * Compute the maximum wave speed of each equation locally.
                                 */
                                
                                std::vector<double> alpha_y;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    alpha_y.push_back(fmax(fabs(lambda_y[ei][idx_cell_B]),
                                        fabs(lambda_y[ei][idx_cell_T])));
                                }
                                
                                boost::multi_array<double, 2> G_y_plus_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                boost::multi_array<double, 2> G_y_minus_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                /*
                                 * Project fluxes onto characteristic fields.
                                 */
                                
                                boost::multi_array<const double*, 2> R_y_inv_intercell(
                                    boost::extents[d_num_eqn][d_num_eqn]);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    for (int ej = 0; ej < d_num_eqn; ej++)
                                    {
                                        R_y_inv_intercell[ei][ej] =
                                            &(projection_matrix->getPointer(1, ei + ej*d_num_eqn)[idx_face_y]);
                                    }
                                }
                                
                                for (int m = 0; m < 6; m++)
                                {
                                    const int idx_cell = (i + d_num_ghosts[0]) +
                                        (j - 3 + m + d_num_ghosts[1])*ghostcell_dims[0];
                                    
                                    std::vector<const double*> F_y_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        F_y_ptr.push_back(&F_y_node[ei][idx_cell]);
                                    }
                                    
                                    std::vector<double> G_y(d_num_eqn, 0.0);
                                    std::vector<double*> G_y_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        G_y_ptr.push_back(&G_y[ei]);
                                    }
                                    
                                    projectVariablesToCharacteristicFields(
                                        G_y_ptr,
                                        F_y_ptr,
                                        R_y_inv_intercell);
                                    
                                    std::vector<const double*> Q_ptr;
                                    Q_ptr.push_back(&rho[idx_cell]);
                                    Q_ptr.push_back(&rho_u[idx_cell]);
                                    Q_ptr.push_back(&rho_v[idx_cell]);
                                    Q_ptr.push_back(&E[idx_cell]);
                                    
                                    std::vector<double> V(d_num_eqn, 0.0);
                                    std::vector<double*> V_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        V_ptr.push_back(&V[ei]);
                                    }
                                    
                                    projectVariablesToCharacteristicFields(
                                        V_ptr,
                                        Q_ptr,
                                        R_y_inv_intercell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        G_y_plus_array[m][ei] = 0.5*(G_y[ei] +
                                            alpha_y[ei]*V[ei]);
                                        
                                        G_y_minus_array[m][ei] = 0.5*(G_y[ei] -
                                            alpha_y[ei]*V[ei]);
                                    }
                                }
                                
                                /*
                                 * Do WENO reconstruction on characteristic variables to get G_y_plus
                                 * and G_y_minus
                                 */
                                
                                std::vector<double> G_y_plus;
                                std::vector<double> G_y_minus;
                                
                                performWENOReconstruction(
                                    G_y_plus,
                                    G_y_minus,
                                    G_y_plus_array,
                                    G_y_minus_array,
                                    Y_DIRECTION);
                                
                                /*
                                 * Project characteristic variables back to physical fields.
                                 */
                                
                                boost::multi_array<const double*, 2> R_y_intercell(
                                    boost::extents[d_num_eqn][d_num_eqn]);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    for (int ej = 0; ej < d_num_eqn; ej++)
                                    {
                                        R_y_intercell[ei][ej] =
                                            &(projection_matrix_inv->getPointer(1, ei + ej*d_num_eqn)[idx_face_y]);
                                    }
                                }
                                
                                std::vector<double> F_y_plus(d_num_eqn, 0.0);
                                std::vector<double> F_y_minus(d_num_eqn, 0.0);
                                
                                std::vector<double*> F_y_plus_ptr;
                                std::vector<double*> F_y_minus_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_y_plus_ptr.push_back(&F_y_plus[ei]);
                                    F_y_minus_ptr.push_back(&F_y_minus[ei]);
                                }
                                
                                std::vector<const double*> G_y_plus_ptr;
                                std::vector<const double*> G_y_minus_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    G_y_plus_ptr.push_back(&G_y_plus[ei]);
                                    G_y_minus_ptr.push_back(&G_y_minus[ei]);
                                }
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    F_y_plus_ptr,
                                    G_y_plus_ptr,
                                    R_y_intercell);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    F_y_minus_ptr,
                                    G_y_minus_ptr,
                                    R_y_intercell);
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double* F_y = convective_flux->getPointer(1, ei);
                                    F_y[idx_face_y] = dt*(F_y_plus[ei] + F_y_minus[ei]);
                                }
                            }
                        }
                    }
                    else if (d_dim == tbox::Dimension(3))
                    {
                        // Get the arrays of time-dependent variables.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* rho_v = momentum->getPointer(1);
                        double* rho_w = momentum->getPointer(2);
                        double* E     = total_energy->getPointer(0);
                        
                        // Get the arrays of temporary patch data.
                        double* u     = velocity->getPointer(0);
                        double* v     = velocity->getPointer(1);
                        double* w     = velocity->getPointer(2);
                        double* p     = pressure->getPointer(0);
                        double* c     = sound_speed->getPointer(0);
                        double* H     = total_specific_enthalpy->getPointer(0);
                        std::vector<double*> F_x_node;
                        std::vector<double*> F_y_node;
                        std::vector<double*> F_z_node;
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
                            F_y_node.push_back(convective_flux_node[1]->getPointer(ei));
                            F_z_node.push_back(convective_flux_node[2]->getPointer(ei));
                        }
                        std::vector<double*> lambda_x;
                        std::vector<double*> lambda_y;
                        std::vector<double*> lambda_z;
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            lambda_x.push_back(wave_speed[0]->getPointer(ei));
                            lambda_y.push_back(wave_speed[1]->getPointer(ei));
                            lambda_z.push_back(wave_speed[2]->getPointer(ei));
                        }
                        
                        // Compute the field of velocities, pressure, sound speed, total specific enthalpy, wave speeds
                        // and fluxes.
                        for (int k = -d_num_ghosts[2]; k < interior_dims[2] + d_num_ghosts[2]; k++)
                        {
                            for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                            {
                                for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                                {
                                    // Compute index into linear data array.
                                    const int idx = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    u[idx] = rho_u[idx]/rho[idx];
                                    v[idx] = rho_v[idx]/rho[idx];
                                    w[idx] = rho_w[idx]/rho[idx];
                                    
                                    std::vector<const double*> m_ptr;
                                    m_ptr.push_back(&rho_u[idx]);
                                    m_ptr.push_back(&rho_v[idx]);
                                    m_ptr.push_back(&rho_w[idx]);
                                    
                                    p[idx] = d_equation_of_state->
                                        getPressure(
                                            &rho[idx],
                                            m_ptr,
                                            &E[idx]);
                                    
                                    c[idx] = d_equation_of_state->
                                        getSoundSpeedWithPressure(
                                            &rho[idx],
                                            &p[idx]);
                                    
                                    H[idx] = (E[idx] + p[idx])/rho[idx];
                                    
                                    F_x_node[0][idx] = rho_u[idx];
                                    F_x_node[1][idx] = rho_u[idx]*u[idx] + p[idx];
                                    F_x_node[2][idx] = rho_u[idx]*v[idx];
                                    F_x_node[3][idx] = rho_u[idx]*w[idx];
                                    F_x_node[4][idx] = u[idx]*(E[idx] + p[idx]);
                                    
                                    F_y_node[0][idx] = rho_v[idx];
                                    F_y_node[1][idx] = rho_v[idx]*u[idx];
                                    F_y_node[2][idx] = rho_v[idx]*v[idx] + p[idx];
                                    F_y_node[3][idx] = rho_v[idx]*w[idx];
                                    F_y_node[4][idx] = v[idx]*(E[idx] + p[idx]);
                                    
                                    F_z_node[0][idx] = rho_w[idx];
                                    F_z_node[1][idx] = rho_w[idx]*u[idx];
                                    F_z_node[2][idx] = rho_w[idx]*v[idx];
                                    F_z_node[3][idx] = rho_w[idx]*w[idx] + p[idx];
                                    F_z_node[4][idx] = w[idx]*(E[idx] + p[idx]);
                                    
                                    lambda_x[0][idx] = u[idx] - c[idx];
                                    lambda_x[1][idx] = u[idx];
                                    lambda_x[2][idx] = u[idx];
                                    lambda_x[3][idx] = u[idx];
                                    lambda_x[4][idx] = u[idx] + c[idx];
                                    
                                    lambda_y[0][idx] = v[idx] - c[idx];
                                    lambda_y[1][idx] = v[idx];
                                    lambda_y[2][idx] = v[idx];
                                    lambda_y[3][idx] = v[idx];
                                    lambda_y[4][idx] = v[idx] + c[idx];
                                    
                                    lambda_z[0][idx] = w[idx] - c[idx];
                                    lambda_z[1][idx] = w[idx];
                                    lambda_z[2][idx] = w[idx];
                                    lambda_z[3][idx] = w[idx];
                                    lambda_z[4][idx] = w[idx] + c[idx];
                                }
                            }
                        }
                        
                        /*
                         * Compute the projection matrix and its inverse at the face normal to the
                         * x direction.
                         */
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = 0; j < interior_dims[1]; j++)
                            {
                                for (int i = 0; i < interior_dims[0] + 1; i++)
                                {
                                    // Compute the indices of left cell, right cell and face of
                                    // projection matrix.
                                    const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_R = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_face_x = i +
                                        j*(interior_dims[0] + 1) +
                                        k*(interior_dims[0] + 1)*interior_dims[1];
                                    
                                    // Get left and right quantities.
                                    const double& u_L = u[idx_cell_L];
                                    const double& u_R = u[idx_cell_R];
                                    
                                    const double& v_L = v[idx_cell_L];
                                    const double& v_R = v[idx_cell_R];
                                    
                                    const double& w_L = w[idx_cell_L];
                                    const double& w_R = w[idx_cell_R];
                                    
                                    const double& c_L = c[idx_cell_L];
                                    const double& c_R = c[idx_cell_R];
                                    
                                    const double& H_L = H[idx_cell_L];
                                    const double& H_R = H[idx_cell_R];
                                    
                                    // Compute simply-averaged quantities.
                                    const double u_average = 0.5*(u_L + u_R);
                                    const double v_average = 0.5*(v_L + v_R);
                                    const double w_average = 0.5*(w_L + w_R);
                                    const double c_average = 0.5*(c_L + c_R);
                                    const double H_average = 0.5*(H_L + H_R);
                                    
                                    boost::multi_array<double*, 2> R_x_intercell(
                                        boost::extents[d_num_eqn][d_num_eqn]);
                                    
                                    boost::multi_array<double*, 2> R_x_inv_intercell(
                                        boost::extents[d_num_eqn][d_num_eqn]);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        for (int ej = 0; ej < d_num_eqn; ej++)
                                        {
                                            R_x_intercell[ei][ej] =
                                                &(projection_matrix_inv->getPointer(0, ei + ej*d_num_eqn)[idx_face_x]);
                                            
                                            R_x_inv_intercell[ei][ej] =
                                                &(projection_matrix->getPointer(0, ei + ej*d_num_eqn)[idx_face_x]);
                                        }
                                    }
                                    
                                    *R_x_intercell[0][0] = 1.0;
                                    *R_x_intercell[0][1] = 0.0;
                                    *R_x_intercell[0][2] = 0.0;
                                    *R_x_intercell[0][3] = 1.0;
                                    *R_x_intercell[0][4] = 1.0;
                                    *R_x_intercell[1][0] = u_average - c_average;
                                    *R_x_intercell[1][1] = 0.0;
                                    *R_x_intercell[1][2] = 0.0;
                                    *R_x_intercell[1][3] = u_average;
                                    *R_x_intercell[1][4] = u_average + c_average;
                                    *R_x_intercell[2][0] = v_average;
                                    *R_x_intercell[2][1] = 1.0;
                                    *R_x_intercell[2][2] = 0.0;
                                    *R_x_intercell[2][3] = v_average;
                                    *R_x_intercell[2][4] = v_average;
                                    *R_x_intercell[3][0] = w_average;
                                    *R_x_intercell[3][1] = 0.0;
                                    *R_x_intercell[3][2] = 1.0;
                                    *R_x_intercell[3][3] = w_average;
                                    *R_x_intercell[3][4] = w_average;
                                    *R_x_intercell[4][0] = H_average - u_average*c_average;
                                    *R_x_intercell[4][1] = v_average;
                                    *R_x_intercell[4][2] = w_average;
                                    *R_x_intercell[4][3] = 0.5*(u_average*u_average + v_average*v_average +
                                        w_average*w_average);
                                    *R_x_intercell[4][4] = H_average + u_average*c_average;
                                    
                                    const double gamma = d_equation_of_state->
                                        getSpeciesThermodynamicProperty(
                                            "gamma",
                                            0);
                                    const double b_1 = (gamma - 1.0)/(c_average*c_average);
                                    const double b_2 = 0.5*(u_average*u_average + v_average*v_average +
                                        w_average*w_average)*b_1;
                                    
                                    *R_x_inv_intercell[0][0] = 0.5*(b_2 + u_average/c_average);
                                    *R_x_inv_intercell[0][1] = -0.5*(b_1*u_average + 1.0/c_average);
                                    *R_x_inv_intercell[0][2] = -0.5*b_1*v_average;
                                    *R_x_inv_intercell[0][3] = -0.5*b_1*w_average;
                                    *R_x_inv_intercell[0][4] = 0.5*b_1;
                                    *R_x_inv_intercell[1][0] = -v_average;
                                    *R_x_inv_intercell[1][1] = 0.0;
                                    *R_x_inv_intercell[1][2] = 1.0;
                                    *R_x_inv_intercell[1][3] = 0.0;
                                    *R_x_inv_intercell[1][4] = 0.0;
                                    *R_x_inv_intercell[2][0] = -w_average;
                                    *R_x_inv_intercell[2][1] = 0.0;
                                    *R_x_inv_intercell[2][2] = 0.0;
                                    *R_x_inv_intercell[2][3] = 1.0;
                                    *R_x_inv_intercell[2][4] = 0.0;
                                    *R_x_inv_intercell[3][0] = 1 - b_2;
                                    *R_x_inv_intercell[3][1] = b_1*u_average;
                                    *R_x_inv_intercell[3][2] = b_1*v_average;
                                    *R_x_inv_intercell[3][3] = b_1*w_average;
                                    *R_x_inv_intercell[3][4] = -b_1;
                                    *R_x_inv_intercell[4][0] = 0.5*(b_2 - u_average/c_average);
                                    *R_x_inv_intercell[4][1] = 0.5*(-b_1*u_average + 1.0/c_average);
                                    *R_x_inv_intercell[4][2] = -0.5*b_1*v_average;
                                    *R_x_inv_intercell[4][3] = -0.5*b_1*w_average;
                                    *R_x_inv_intercell[4][4] = 0.5*b_1;
                                }
                            }
                        }
                        
                        /*
                         * Compute the projection matrix and its inverse at the face normal to the
                         * y direction.
                         */
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int k = 0; k < interior_dims[2]; k++)
                            {
                                for (int j = 0; j < interior_dims[1] + 1; j++)
                                {
                                    // Compute the indices of bottom cell, top cell and face of
                                    // projection matrix.
                                    const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_T = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_face_y = j +
                                        k*(interior_dims[1] + 1) +
                                        i*(interior_dims[1] + 1)*interior_dims[2];
                                    
                                    // Get bottom and top quantities.
                                    const double& u_B = u[idx_cell_B];
                                    const double& u_T = u[idx_cell_T];
                                    
                                    const double& v_B = v[idx_cell_B];
                                    const double& v_T = v[idx_cell_T];
                                    
                                    const double& w_B = w[idx_cell_B];
                                    const double& w_T = w[idx_cell_T];
                                    
                                    const double& c_B = c[idx_cell_B];
                                    const double& c_T = c[idx_cell_T];
                                    
                                    const double& H_B = H[idx_cell_B];
                                    const double& H_T = H[idx_cell_T];
                                    
                                    // Compute simply-averaged quantities.
                                    const double u_average = 0.5*(u_B + u_T);
                                    const double v_average = 0.5*(v_B + v_T);
                                    const double w_average = 0.5*(w_B + w_T);
                                    const double c_average = 0.5*(c_B + c_T);
                                    const double H_average = 0.5*(H_B + H_T);
                                    
                                    boost::multi_array<double*, 2> R_y_intercell(
                                        boost::extents[d_num_eqn][d_num_eqn]);
                                    
                                    boost::multi_array<double*, 2> R_y_inv_intercell(
                                        boost::extents[d_num_eqn][d_num_eqn]);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        for (int ej = 0; ej < d_num_eqn; ej++)
                                        {
                                            R_y_intercell[ei][ej] =
                                                &(projection_matrix_inv->getPointer(1, ei + ej*d_num_eqn)[idx_face_y]);
                                            
                                            R_y_inv_intercell[ei][ej] =
                                                &(projection_matrix->getPointer(1, ei + ej*d_num_eqn)[idx_face_y]);
                                        }
                                    }
                                    
                                    *R_y_intercell[0][0] = 1.0;
                                    *R_y_intercell[0][1] = 0.0;
                                    *R_y_intercell[0][2] = 0.0;
                                    *R_y_intercell[0][3] = 1.0;
                                    *R_y_intercell[0][4] = 1.0;
                                    *R_y_intercell[1][0] = u_average;
                                    *R_y_intercell[1][1] = 1.0;
                                    *R_y_intercell[1][2] = 0.0;
                                    *R_y_intercell[1][3] = u_average;
                                    *R_y_intercell[1][4] = u_average;
                                    *R_y_intercell[2][0] = v_average - c_average;
                                    *R_y_intercell[2][1] = 0.0;
                                    *R_y_intercell[2][2] = 0.0;
                                    *R_y_intercell[2][3] = v_average;
                                    *R_y_intercell[2][4] = v_average + c_average;
                                    *R_y_intercell[3][0] = w_average;
                                    *R_y_intercell[3][1] = 0.0;
                                    *R_y_intercell[3][2] = 1.0;
                                    *R_y_intercell[3][3] = w_average;
                                    *R_y_intercell[3][4] = w_average;
                                    *R_y_intercell[4][0] = H_average - v_average*c_average;
                                    *R_y_intercell[4][1] = u_average;
                                    *R_y_intercell[4][2] = w_average;
                                    *R_y_intercell[4][3] = 0.5*(u_average*u_average + v_average*v_average +
                                        w_average*w_average);
                                    *R_y_intercell[4][4] = H_average + v_average*c_average;
                                    
                                    const double gamma = d_equation_of_state->
                                        getSpeciesThermodynamicProperty(
                                            "gamma",
                                            0);
                                    const double b_1 = (gamma - 1)/(c_average*c_average);
                                    const double b_2 = 0.5*(u_average*u_average + v_average*v_average +
                                        w_average*w_average)*b_1;
                                    
                                    *R_y_inv_intercell[0][0] = 0.5*(b_2 + v_average/c_average);
                                    *R_y_inv_intercell[0][1] = -0.5*b_1*u_average;
                                    *R_y_inv_intercell[0][2] = -0.5*(b_1*v_average + 1.0/c_average);
                                    *R_y_inv_intercell[0][3] = -0.5*b_1*w_average;
                                    *R_y_inv_intercell[0][4] = 0.5*b_1;
                                    *R_y_inv_intercell[1][0] = -u_average;
                                    *R_y_inv_intercell[1][1] = 1.0;
                                    *R_y_inv_intercell[1][2] = 0.0;
                                    *R_y_inv_intercell[1][3] = 0.0;
                                    *R_y_inv_intercell[1][4] = 0.0;
                                    *R_y_inv_intercell[2][0] = -w_average;
                                    *R_y_inv_intercell[2][1] = 0.0;
                                    *R_y_inv_intercell[2][2] = 0.0;
                                    *R_y_inv_intercell[2][3] = 1.0;
                                    *R_y_inv_intercell[2][4] = 0.0;
                                    *R_y_inv_intercell[3][0] = 1 - b_2;
                                    *R_y_inv_intercell[3][1] = b_1*u_average;
                                    *R_y_inv_intercell[3][2] = b_1*v_average;
                                    *R_y_inv_intercell[3][3] = b_1*w_average;
                                    *R_y_inv_intercell[3][4] = -b_1;
                                    *R_y_inv_intercell[4][0] = 0.5*(b_2 - v_average/c_average);
                                    *R_y_inv_intercell[4][1] = -0.5*b_1*u_average;
                                    *R_y_inv_intercell[4][2] = 0.5*(-b_1*v_average + 1.0/c_average);
                                    *R_y_inv_intercell[4][3] = -0.5*b_1*w_average;
                                    *R_y_inv_intercell[4][4] = 0.5*b_1;
                                }
                            }
                        }
                        
                        /*
                         * Compute the projection matrix and its inverse at the face normal to the
                         * z direction.
                         */
                        // Compute the fluxes in the z direction.
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0]; i++)
                            {
                                for (int k = 0; k < interior_dims[2] + 1; k++)
                                {
                                    // Compute the indices of bottom cell, top cell and face of
                                    // projection matrix.
                                    const int idx_cell_B = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_F = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_face_z = k +
                                        i*(interior_dims[2] + 1) +
                                        j*(interior_dims[2] + 1)*interior_dims[0];
                                    
                                    // Get back and front quantities.
                                    const double& u_B = u[idx_cell_B];
                                    const double& u_F = u[idx_cell_F];
                                    
                                    const double& v_B = v[idx_cell_B];
                                    const double& v_F = v[idx_cell_F];
                                    
                                    const double& w_B = w[idx_cell_B];
                                    const double& w_F = w[idx_cell_F];
                                    
                                    const double& c_B = c[idx_cell_B];
                                    const double& c_F = c[idx_cell_F];
                                    
                                    const double& H_B = H[idx_cell_B];
                                    const double& H_F = H[idx_cell_F];
                                    
                                    // Compute simply-averaged quantities.
                                    const double u_average = 0.5*(u_B + u_F);
                                    const double v_average = 0.5*(v_B + v_F);
                                    const double w_average = 0.5*(w_B + w_F);
                                    const double c_average = 0.5*(c_B + c_F);
                                    const double H_average = 0.5*(H_B + H_F);
                                    
                                    boost::multi_array<double*, 2> R_z_intercell(
                                        boost::extents[d_num_eqn][d_num_eqn]);
                                    
                                    boost::multi_array<double*, 2> R_z_inv_intercell(
                                        boost::extents[d_num_eqn][d_num_eqn]);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        for (int ej = 0; ej < d_num_eqn; ej++)
                                        {
                                            R_z_intercell[ei][ej] =
                                                &(projection_matrix_inv->getPointer(2, ei + ej*d_num_eqn)[idx_face_z]);
                                            
                                            R_z_inv_intercell[ei][ej] =
                                                &(projection_matrix->getPointer(2, ei + ej*d_num_eqn)[idx_face_z]);
                                        }
                                    }
                                    
                                    *R_z_intercell[0][0] = 1.0;
                                    *R_z_intercell[0][1] = 0.0;
                                    *R_z_intercell[0][2] = 0.0;
                                    *R_z_intercell[0][3] = 1.0;
                                    *R_z_intercell[0][4] = 1.0;
                                    *R_z_intercell[1][0] = u_average;
                                    *R_z_intercell[1][1] = 1.0;
                                    *R_z_intercell[1][2] = 0.0;
                                    *R_z_intercell[1][3] = u_average;
                                    *R_z_intercell[1][4] = u_average;
                                    *R_z_intercell[2][0] = v_average;
                                    *R_z_intercell[2][1] = 0.0;
                                    *R_z_intercell[2][2] = 1.0;
                                    *R_z_intercell[2][3] = v_average;
                                    *R_z_intercell[2][4] = v_average;
                                    *R_z_intercell[3][0] = w_average - c_average;
                                    *R_z_intercell[3][1] = 0.0;
                                    *R_z_intercell[3][2] = 0.0;
                                    *R_z_intercell[3][3] = w_average;
                                    *R_z_intercell[3][4] = w_average + c_average;
                                    *R_z_intercell[4][0] = H_average - w_average*c_average;
                                    *R_z_intercell[4][1] = u_average;
                                    *R_z_intercell[4][2] = v_average;
                                    *R_z_intercell[4][3] = 0.5*(u_average*u_average + v_average*v_average +
                                        w_average*w_average);
                                    *R_z_intercell[4][4] = H_average + w_average*c_average;
                                    
                                    const double gamma = d_equation_of_state->
                                        getSpeciesThermodynamicProperty(
                                            "gamma",
                                            0);
                                    const double b_1 = (gamma - 1)/(c_average*c_average);
                                    const double b_2 = 0.5*(u_average*u_average + v_average*v_average +
                                        w_average*w_average)*b_1;
                                    
                                    *R_z_inv_intercell[0][0] = 0.5*(b_2 + w_average/c_average);
                                    *R_z_inv_intercell[0][1] = -0.5*b_1*u_average;
                                    *R_z_inv_intercell[0][2] = -0.5*b_1*v_average;
                                    *R_z_inv_intercell[0][3] = -0.5*(b_1*w_average + 1.0/c_average);
                                    *R_z_inv_intercell[0][4] = 0.5*b_1;
                                    *R_z_inv_intercell[1][0] = -u_average;
                                    *R_z_inv_intercell[1][1] = 1.0;
                                    *R_z_inv_intercell[1][2] = 0.0;
                                    *R_z_inv_intercell[1][3] = 0.0;
                                    *R_z_inv_intercell[1][4] = 0.0;
                                    *R_z_inv_intercell[2][0] = -v_average;
                                    *R_z_inv_intercell[2][1] = 0.0;
                                    *R_z_inv_intercell[2][2] = 1.0;
                                    *R_z_inv_intercell[2][3] = 0.0;
                                    *R_z_inv_intercell[2][4] = 0.0;
                                    *R_z_inv_intercell[3][0] = 1 - b_2;
                                    *R_z_inv_intercell[3][1] = b_1*u_average;
                                    *R_z_inv_intercell[3][2] = b_1*v_average;
                                    *R_z_inv_intercell[3][3] = b_1*w_average;
                                    *R_z_inv_intercell[3][4] = -b_1;
                                    *R_z_inv_intercell[4][0] = 0.5*(b_2 - w_average/c_average);
                                    *R_z_inv_intercell[4][1] = -0.5*b_1*u_average;
                                    *R_z_inv_intercell[4][2] = -0.5*b_1*v_average;
                                    *R_z_inv_intercell[4][3] = -0.5*(b_1*w_average - 1.0/c_average);
                                    *R_z_inv_intercell[4][4] = 0.5*b_1;
                                }
                            }
                        }
                        
                        // Compute the fluxes in the x direction.
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = 0; j < interior_dims[1]; j++)
                            {
                                for (int i = 0; i < interior_dims[0] + 1; i++)
                                {
                                    // Compute the indices of left cell, right cell and face.
                                    const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_R = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_face_x = i +
                                        j*(interior_dims[0] + 1) +
                                        k*(interior_dims[0] + 1)*interior_dims[1];
                                    
                                    /*
                                     * Compute the maximum wave speed of each equation locally.
                                     */
                                    
                                    std::vector<double> alpha_x;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        alpha_x.push_back(fmax(fabs(lambda_x[ei][idx_cell_L]),
                                            fabs(lambda_x[ei][idx_cell_R])));
                                    }
                                    
                                    boost::multi_array<double, 2> G_x_plus_array(
                                        boost::extents[6][d_num_eqn]);
                                    
                                    boost::multi_array<double, 2> G_x_minus_array(
                                        boost::extents[6][d_num_eqn]);
                                    
                                    /*
                                     * Project fluxes onto characteristic fields.
                                     */
                                    
                                    boost::multi_array<const double*, 2> R_x_inv_intercell(
                                        boost::extents[d_num_eqn][d_num_eqn]);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        for (int ej = 0; ej < d_num_eqn; ej++)
                                        {
                                            R_x_inv_intercell[ei][ej] =
                                                &(projection_matrix->getPointer(0, ei + ej*d_num_eqn)[idx_face_x]);
                                        }
                                    }
                                    
                                    for (int m = 0; m < 6; m++)
                                    {
                                        const int idx_cell = (i - 3 + m + d_num_ghosts[0]) +
                                            (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                            (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                        
                                        std::vector<const double*> F_x_ptr;
                                        for (int ei = 0; ei < d_num_eqn; ei++)
                                        {
                                            F_x_ptr.push_back(&F_x_node[ei][idx_cell]);
                                        }
                                        
                                        std::vector<double> G_x(d_num_eqn, 0.0);
                                        std::vector<double*> G_x_ptr;
                                        for (int ei = 0; ei < d_num_eqn; ei++)
                                        {
                                            G_x_ptr.push_back(&G_x[ei]);
                                        }
                                        
                                        projectVariablesToCharacteristicFields(
                                            G_x_ptr,
                                            F_x_ptr,
                                            R_x_inv_intercell);
                                        
                                        std::vector<const double*> Q_ptr;
                                        Q_ptr.push_back(&rho[idx_cell]);
                                        Q_ptr.push_back(&rho_u[idx_cell]);
                                        Q_ptr.push_back(&rho_v[idx_cell]);
                                        Q_ptr.push_back(&rho_w[idx_cell]);
                                        Q_ptr.push_back(&E[idx_cell]);
                                        
                                        std::vector<double> V(d_num_eqn, 0.0);
                                        std::vector<double*> V_ptr;
                                        for (int ei = 0; ei < d_num_eqn; ei++)
                                        {
                                            V_ptr.push_back(&V[ei]);
                                        }
                                        
                                        projectVariablesToCharacteristicFields(
                                            V_ptr,
                                            Q_ptr,
                                            R_x_inv_intercell);
                                        
                                        for (int ei = 0; ei < d_num_eqn; ei++)
                                        {
                                            G_x_plus_array[m][ei] = 0.5*(G_x[ei] +
                                                alpha_x[ei]*V[ei]);
                                            
                                            G_x_minus_array[m][ei] = 0.5*(G_x[ei] -
                                                alpha_x[ei]*V[ei]);
                                        }
                                    }
                                    
                                    /*
                                     * Do WENO reconstruction on characteristic variables to get G_x_plus
                                     * and G_x_minus
                                     */
                                    
                                    std::vector<double> G_x_plus;
                                    std::vector<double> G_x_minus;
                                    
                                    performWENOReconstruction(
                                        G_x_plus,
                                        G_x_minus,
                                        G_x_plus_array,
                                        G_x_minus_array,
                                        X_DIRECTION);
                                    
                                    /*
                                     * Project characteristic variables back to physical fields.
                                     */
                                    
                                    boost::multi_array<const double*, 2> R_x_intercell(
                                        boost::extents[d_num_eqn][d_num_eqn]);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        for (int ej = 0; ej < d_num_eqn; ej++)
                                        {
                                            R_x_intercell[ei][ej] =
                                                &(projection_matrix_inv->getPointer(0, ei + ej*d_num_eqn)[idx_face_x]);
                                        }
                                    }
                                    
                                    std::vector<double> F_x_plus(d_num_eqn, 0.0);
                                    std::vector<double> F_x_minus(d_num_eqn, 0.0);
                                    
                                    std::vector<double*> F_x_plus_ptr;
                                    std::vector<double*> F_x_minus_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        F_x_plus_ptr.push_back(&F_x_plus[ei]);
                                        F_x_minus_ptr.push_back(&F_x_minus[ei]);
                                    }
                                    
                                    std::vector<const double*> G_x_plus_ptr;
                                    std::vector<const double*> G_x_minus_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        G_x_plus_ptr.push_back(&G_x_plus[ei]);
                                        G_x_minus_ptr.push_back(&G_x_minus[ei]);
                                    }
                                    
                                    projectCharacteristicVariablesToPhysicalFields(
                                        F_x_plus_ptr,
                                        G_x_plus_ptr,
                                        R_x_intercell);
                                    
                                    projectCharacteristicVariablesToPhysicalFields(
                                        F_x_minus_ptr,
                                        G_x_minus_ptr,
                                        R_x_intercell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_x = convective_flux->getPointer(0, ei);
                                        F_x[idx_face_x] = dt*(F_x_plus[ei] + F_x_minus[ei]);
                                    }
                                }
                            }
                        }
                        
                        // Compute the fluxes in the y direction.
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int k = 0; k < interior_dims[2]; k++)
                            {
                                for (int j = 0; j < interior_dims[1] + 1; j++)
                                {
                                    // Compute the indices of bottom cell, top cell and face.
                                    const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_T = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_face_y = j +
                                        k*(interior_dims[1] + 1) +
                                        i*(interior_dims[1] + 1)*interior_dims[2];
                                    
                                    /*
                                     * Compute the maximum wave speed of each equation locally.
                                     */
                                    
                                    std::vector<double> alpha_y;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        alpha_y.push_back(fmax(fabs(lambda_y[ei][idx_cell_B]),
                                            fabs(lambda_y[ei][idx_cell_T])));
                                    }
                                    
                                    boost::multi_array<double, 2> G_y_plus_array(
                                        boost::extents[6][d_num_eqn]);
                                    
                                    boost::multi_array<double, 2> G_y_minus_array(
                                        boost::extents[6][d_num_eqn]);
                                    
                                    /*
                                     * Project fluxes onto characteristic fields.
                                     */
                                    
                                    boost::multi_array<const double*, 2> R_y_inv_intercell(
                                        boost::extents[d_num_eqn][d_num_eqn]);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        for (int ej = 0; ej < d_num_eqn; ej++)
                                        {
                                            R_y_inv_intercell[ei][ej] =
                                                &(projection_matrix->getPointer(1, ei + ej*d_num_eqn)[idx_face_y]);
                                        }
                                    }
                                    
                                    for (int m = 0; m < 6; m++)
                                    {
                                        const int idx_cell = (i + d_num_ghosts[0]) +
                                            (j - 3 + m + d_num_ghosts[1])*ghostcell_dims[0] +
                                            (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                        
                                        std::vector<const double*> F_y_ptr;
                                        for (int ei = 0; ei < d_num_eqn; ei++)
                                        {
                                            F_y_ptr.push_back(&F_y_node[ei][idx_cell]);
                                        }
                                        
                                        std::vector<double> G_y(d_num_eqn, 0.0);
                                        std::vector<double*> G_y_ptr;
                                        for (int ei = 0; ei < d_num_eqn; ei++)
                                        {
                                            G_y_ptr.push_back(&G_y[ei]);
                                        }
                                        
                                        projectVariablesToCharacteristicFields(
                                            G_y_ptr,
                                            F_y_ptr,
                                            R_y_inv_intercell);
                                        
                                        std::vector<const double*> Q_ptr;
                                        Q_ptr.push_back(&rho[idx_cell]);
                                        Q_ptr.push_back(&rho_u[idx_cell]);
                                        Q_ptr.push_back(&rho_v[idx_cell]);
                                        Q_ptr.push_back(&rho_w[idx_cell]);
                                        Q_ptr.push_back(&E[idx_cell]);
                                        
                                        std::vector<double> V(d_num_eqn, 0.0);
                                        std::vector<double*> V_ptr;
                                        for (int ei = 0; ei < d_num_eqn; ei++)
                                        {
                                            V_ptr.push_back(&V[ei]);
                                        }
                                        
                                        projectVariablesToCharacteristicFields(
                                            V_ptr,
                                            Q_ptr,
                                            R_y_inv_intercell);
                                        
                                        for (int ei = 0; ei < d_num_eqn; ei++)
                                        {
                                            G_y_plus_array[m][ei] = 0.5*(G_y[ei] +
                                                alpha_y[ei]*V[ei]);
                                            
                                            G_y_minus_array[m][ei] = 0.5*(G_y[ei] -
                                                alpha_y[ei]*V[ei]);
                                        }
                                    }
                                    
                                    /*
                                     * Do WENO reconstruction on characteristic variables to get G_y_plus
                                     * and G_y_minus
                                     */
                                    
                                    std::vector<double> G_y_plus;
                                    std::vector<double> G_y_minus;
                                    
                                    performWENOReconstruction(
                                        G_y_plus, G_y_minus,
                                        G_y_plus_array,
                                        G_y_minus_array,
                                        Y_DIRECTION);
                                    
                                    /*
                                     * Project characteristic variables back to physical fields.
                                     */
                                    
                                    boost::multi_array<const double*, 2> R_y_intercell(
                                        boost::extents[d_num_eqn][d_num_eqn]);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        for (int ej = 0; ej < d_num_eqn; ej++)
                                        {
                                            R_y_intercell[ei][ej] =
                                                &(projection_matrix_inv->getPointer(1, ei + ej*d_num_eqn)[idx_face_y]);
                                        }
                                    }
                                    
                                    std::vector<double> F_y_plus(d_num_eqn, 0.0);
                                    std::vector<double> F_y_minus(d_num_eqn, 0.0);
                                    
                                    std::vector<double*> F_y_plus_ptr;
                                    std::vector<double*> F_y_minus_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        F_y_plus_ptr.push_back(&F_y_plus[ei]);
                                        F_y_minus_ptr.push_back(&F_y_minus[ei]);
                                    }
                                    
                                    std::vector<const double*> G_y_plus_ptr;
                                    std::vector<const double*> G_y_minus_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        G_y_plus_ptr.push_back(&G_y_plus[ei]);
                                        G_y_minus_ptr.push_back(&G_y_minus[ei]);
                                    }
                                    
                                    projectCharacteristicVariablesToPhysicalFields(
                                        F_y_plus_ptr,
                                        G_y_plus_ptr,
                                        R_y_intercell);
                                    
                                    projectCharacteristicVariablesToPhysicalFields(
                                        F_y_minus_ptr,
                                        G_y_minus_ptr,
                                        R_y_intercell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_y = convective_flux->getPointer(1, ei);
                                        F_y[idx_face_y] = dt*(F_y_plus[ei] + F_y_minus[ei]);
                                    }
                                }
                            }
                        }
                        
                        // Compute the fluxes in the z direction.
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0]; i++)
                            {
                                for (int k = 0; k < interior_dims[2] + 1; k++)
                                {
                                    // Compute the indices of back cell, front cell and face.
                                    const int idx_cell_B = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_F = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_face_z = k +
                                        i*(interior_dims[2] + 1) +
                                        j*(interior_dims[2] + 1)*interior_dims[0];
                                    
                                    /*
                                     * Compute the maximum wave speed of each equation locally.
                                     */
                                    
                                    std::vector<double> alpha_z;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        alpha_z.push_back(fmax(fabs(lambda_z[ei][idx_cell_B]),
                                            fabs(lambda_z[ei][idx_cell_F])));
                                    }
                                    
                                    boost::multi_array<double, 2> G_z_plus_array(
                                        boost::extents[6][d_num_eqn]);
                                    
                                    boost::multi_array<double, 2> G_z_minus_array(
                                        boost::extents[6][d_num_eqn]);
                                    
                                    /*
                                     * Project fluxes onto characteristic fields.
                                     */
                                    
                                    boost::multi_array<const double*, 2> R_z_inv_intercell(
                                        boost::extents[d_num_eqn][d_num_eqn]);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        for (int ej = 0; ej < d_num_eqn; ej++)
                                        {
                                            R_z_inv_intercell[ei][ej] =
                                                &(projection_matrix->getPointer(2, ei + ej*d_num_eqn)[idx_face_z]);
                                        }
                                    }
                                    
                                    for (int m = 0; m < 6; m++)
                                    {
                                        const int idx_cell = (i + d_num_ghosts[0]) +
                                            (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                            (k - 3 + m + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                        
                                        std::vector<const double*> F_z_ptr;
                                        for (int ei = 0; ei < d_num_eqn; ei++)
                                        {
                                            F_z_ptr.push_back(&F_z_node[ei][idx_cell]);
                                        }
                                        
                                        std::vector<double> G_z(d_num_eqn, 0.0);
                                        std::vector<double*> G_z_ptr;
                                        for (int ei = 0; ei < d_num_eqn; ei++)
                                        {
                                            G_z_ptr.push_back(&G_z[ei]);
                                        }
                                        
                                        projectVariablesToCharacteristicFields(
                                            G_z_ptr,
                                            F_z_ptr,
                                            R_z_inv_intercell);
                                        
                                        std::vector<const double*> Q_ptr;
                                        Q_ptr.push_back(&rho[idx_cell]);
                                        Q_ptr.push_back(&rho_u[idx_cell]);
                                        Q_ptr.push_back(&rho_v[idx_cell]);
                                        Q_ptr.push_back(&rho_w[idx_cell]);
                                        Q_ptr.push_back(&E[idx_cell]);
                                        
                                        std::vector<double> V(d_num_eqn, 0.0);
                                        std::vector<double*> V_ptr;
                                        for (int ei = 0; ei < d_num_eqn; ei++)
                                        {
                                            V_ptr.push_back(&V[ei]);
                                        }
                                        
                                        projectVariablesToCharacteristicFields(
                                            V_ptr,
                                            Q_ptr,
                                            R_z_inv_intercell);
                                        
                                        for (int ei = 0; ei < d_num_eqn; ei++)
                                        {
                                            G_z_plus_array[m][ei] = 0.5*(G_z[ei] +
                                                alpha_z[ei]*V[ei]);
                                            
                                            G_z_minus_array[m][ei] = 0.5*(G_z[ei] -
                                                alpha_z[ei]*V[ei]);
                                        }
                                    }
                                    
                                    /*
                                     * Do WENO reconstruction on characteristic variables to get G_z_plus
                                     * and G_z_minus
                                     */
                                    
                                    std::vector<double> G_z_plus;
                                    std::vector<double> G_z_minus;
                                    
                                    performWENOReconstruction(
                                        G_z_plus,
                                        G_z_minus,
                                        G_z_plus_array,
                                        G_z_minus_array,
                                        Z_DIRECTION);
                                    
                                    /*
                                     * Project characteristic variables back to physical fields.
                                     */
                                    
                                    boost::multi_array<const double*, 2> R_z_intercell(
                                        boost::extents[d_num_eqn][d_num_eqn]);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        for (int ej = 0; ej < d_num_eqn; ej++)
                                        {
                                            R_z_intercell[ei][ej] =
                                                &(projection_matrix_inv->getPointer(2, ei + ej*d_num_eqn)[idx_face_z]);
                                        }
                                    }
                                    
                                    std::vector<double> F_z_plus(d_num_eqn, 0.0);
                                    std::vector<double> F_z_minus(d_num_eqn, 0.0);
                                    
                                    std::vector<double*> F_z_plus_ptr;
                                    std::vector<double*> F_z_minus_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        F_z_plus_ptr.push_back(&F_z_plus[ei]);
                                        F_z_minus_ptr.push_back(&F_z_minus[ei]);
                                    }
                                    
                                    std::vector<const double*> G_z_plus_ptr;
                                    std::vector<const double*> G_z_minus_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        G_z_plus_ptr.push_back(&G_z_plus[ei]);
                                        G_z_minus_ptr.push_back(&G_z_minus[ei]);
                                    }
                                    
                                    projectCharacteristicVariablesToPhysicalFields(
                                        F_z_plus_ptr,
                                        G_z_plus_ptr,
                                        R_z_intercell);
                                    
                                    projectCharacteristicVariablesToPhysicalFields(
                                        F_z_minus_ptr,
                                        G_z_minus_ptr,
                                        R_z_intercell);
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_z = convective_flux->getPointer(2, ei);
                                        F_z[idx_face_z] = dt*(F_z_plus[ei] + F_z_minus[ei]);
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    TBOX_ERROR(d_object_name
                           << ": "
                           << "The supplied equation of state cannot be used."
                           << " Only ideal gas assumption can be used with"
                           << " the convective flux reconstructor."
                           << std::endl);
                }
                
                break;
            }
            case FOUR_EQN_SHYUE:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "FOUR_EQN_SHYUE flow model is not implemented."
                           << std::endl);
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "FIVE_EQN_ALLAIRE flow model is not implemented."
                           << std::endl);
                
                break;
            }
            default:
            {
                TBOX_ERROR(d_object_name
                           << ": "
                           << "Unknown d_flow_model."
                           << std::endl);
            }
        }
        
    }
    else
    {
        TBOX_ERROR(d_object_name
                   << ": "
                   << "Variables are not set."
                   << std::endl);
    }
}

/*
 * Convert fluxes/conservative variables into characteristic variables.
 */
void
ConvectiveFluxReconstructorWENO_CU6_M2_LLF::projectVariablesToCharacteristicFields(
    std::vector<double*> characteristic_variables,
    const std::vector<const double*> physical_variables,
    const boost::multi_array<const double*, 2> projection_matrix)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(characteristic_variables.size() ==
        physical_variables.size());
    
    TBOX_ASSERT(projection_matrix.num_dimensions() == 2);
    
    TBOX_ASSERT(projection_matrix.shape()[0] ==
        projection_matrix.shape()[1]);
    
    TBOX_ASSERT(projection_matrix.size() ==
        physical_variables.size());
#endif
    
    std::vector<double*>&                         W               = characteristic_variables;
    const std::vector<const double*>&             V               = physical_variables;
    const boost::multi_array<const double*, 2>&   R_inv_intercell = projection_matrix;
    
    for (int ri = 0; ri < static_cast<int>(R_inv_intercell.shape()[0]); ri++)
    {
        *(W[ri]) = 0.0;
        for (int ci = 0; ci < static_cast<int>(R_inv_intercell.shape()[1]); ci++)
        {
            *(W[ri]) += *(R_inv_intercell[ri][ci])*(*(V[ci]));
        }
    }
}


/*
 * Convert characteristic variables into fluxes/conservative variables.
 */
void
ConvectiveFluxReconstructorWENO_CU6_M2_LLF::projectCharacteristicVariablesToPhysicalFields(
   std::vector<double*> physical_variables,
   const std::vector<const double*> characteristic_variables,
   const boost::multi_array<const double*, 2> projection_matrix_inv)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(characteristic_variables.size() ==
        physical_variables.size());
    
    TBOX_ASSERT(projection_matrix_inv.num_dimensions() == 2);
    
    TBOX_ASSERT(projection_matrix_inv.shape()[0] ==
        projection_matrix_inv.shape()[1]);
    
    TBOX_ASSERT(projection_matrix_inv.size() ==
        physical_variables.size());
#endif
    
    std::vector<double*>&                         V           = physical_variables;
    const std::vector<const double*>&             W           = characteristic_variables;
    const boost::multi_array<const double*, 2>&   R_intercell = projection_matrix_inv;
    
    for (int ri = 0; ri < static_cast<int>(R_intercell.shape()[0]); ri++)
    {
        *(V[ri]) = 0.0;
        for (int ci = 0; ci < static_cast<int>(R_intercell.shape()[1]); ci++)
        {
            *(V[ri]) += *(R_intercell[ri][ci])*(*(W[ci]));
        }
    }
}


/*
 * Compute beta_plus's.
 */
boost::multi_array<double, 2>
ConvectiveFluxReconstructorWENO_CU6_M2_LLF::computeBetaPlus(
   const boost::multi_array<double, 2>& G_plus_array)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(G_plus_array.shape()[0]) == 6);
    TBOX_ASSERT(static_cast<int>(G_plus_array.shape()[1]) == d_num_eqn);
#endif
    
    boost::multi_array<double, 2> beta_plus(
        boost::extents[d_num_eqn][4]);
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        beta_plus[ei][0] = 13.0/12*pow((G_plus_array[0][ei] - 2*G_plus_array[1][ei] + G_plus_array[2][ei]), 2) +
            1.0/4*pow((G_plus_array[0][ei] - 4*G_plus_array[1][ei] + 3*G_plus_array[2][ei]), 2);
        
        beta_plus[ei][1] = 13.0/12*pow((G_plus_array[1][ei] - 2*G_plus_array[2][ei] + G_plus_array[3][ei]), 2) +
            1.0/4*pow((G_plus_array[1][ei] - G_plus_array[3][ei]), 2);
        
        beta_plus[ei][2] = 13.0/12*pow((G_plus_array[2][ei] - 2*G_plus_array[3][ei] + G_plus_array[4][ei]), 2) +
            1.0/4*pow((3*G_plus_array[2][ei] - 4*G_plus_array[3][ei] + G_plus_array[4][ei]), 2);
        
        beta_plus[ei][3] = 1.0/120960*(271779*G_plus_array[0][ei]*G_plus_array[0][ei] + G_plus_array[0][ei]*
            (-2380800*G_plus_array[1][ei] + 4086352*G_plus_array[2][ei] - 3462252*G_plus_array[3][ei] +
            1458762*G_plus_array[4][ei] - 245620*G_plus_array[5][ei]) + G_plus_array[1][ei]*
            (5653317*G_plus_array[1][ei] - 20427884*G_plus_array[2][ei] + 17905032*G_plus_array[3][ei] -
            7727988*G_plus_array[4][ei] + 1325006*G_plus_array[5][ei]) + G_plus_array[2][ei]*
            (19510972*G_plus_array[2][ei] - 35817664*G_plus_array[3][ei] + 15929912*G_plus_array[4][ei] -
            2792660*G_plus_array[5][ei]) + G_plus_array[3][ei]*(17195652*G_plus_array[3][ei] -
            15880404*G_plus_array[4][ei] + 2863984*G_plus_array[5][ei]) + G_plus_array[4][ei]*
            (3824847*G_plus_array[4][ei] - 1429976*G_plus_array[5][ei]) +
            139633*G_plus_array[5][ei]*G_plus_array[5][ei]);
    }
    
    return beta_plus;
}


/*
 * Compute beta_minus's.
 */
boost::multi_array<double, 2>
ConvectiveFluxReconstructorWENO_CU6_M2_LLF::computeBetaMinus(
   const boost::multi_array<double, 2>& G_minus_array)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(G_minus_array.shape()[0]) == 6);
    TBOX_ASSERT(static_cast<int>(G_minus_array.shape()[1]) == d_num_eqn);
#endif
    
    boost::multi_array<double, 2> beta_minus(
        boost::extents[d_num_eqn][4]);
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        beta_minus[ei][0] = 13.0/12*pow((G_minus_array[5][ei] - 2*G_minus_array[4][ei] + G_minus_array[3][ei]), 2) +
            1.0/4*pow((G_minus_array[5][ei] - 4*G_minus_array[4][ei] + 3*G_minus_array[3][ei]), 2);
        
        beta_minus[ei][1] = 13.0/12*pow((G_minus_array[4][ei] - 2*G_minus_array[3][ei] + G_minus_array[2][ei]), 2) +
            1.0/4*pow((G_minus_array[4][ei] - G_minus_array[2][ei]), 2);
        
        beta_minus[ei][2] = 13.0/12*pow((G_minus_array[3][ei] - 2*G_minus_array[2][ei] + G_minus_array[1][ei]), 2) +
            1.0/4*pow((3*G_minus_array[3][ei] - 4*G_minus_array[2][ei] + G_minus_array[1][ei]), 2);
        
        beta_minus[ei][3] = 1.0/120960*(271779*G_minus_array[5][ei]*G_minus_array[5][ei] + G_minus_array[5][ei]*
            (-2380800*G_minus_array[4][ei] + 4086352*G_minus_array[3][ei] - 3462252*G_minus_array[2][ei] +
            1458762*G_minus_array[1][ei] - 245620*G_minus_array[0][ei]) + G_minus_array[4][ei]*
            (5653317*G_minus_array[4][ei] - 20427884*G_minus_array[3][ei] + 17905032*G_minus_array[2][ei] -
            7727988*G_minus_array[1][ei] + 1325006*G_minus_array[0][ei]) + G_minus_array[3][ei]*
            (19510972*G_minus_array[3][ei] - 35817664*G_minus_array[2][ei] + 15929912*G_minus_array[1][ei] -
            2792660*G_minus_array[0][ei]) + G_minus_array[2][ei]*(17195652*G_minus_array[2][ei] -
            15880404*G_minus_array[1][ei] + 2863984*G_minus_array[0][ei]) + G_minus_array[1][ei]*
            (3824847*G_minus_array[1][ei] - 1429976*G_minus_array[0][ei]) +
            139633*G_minus_array[0][ei]*G_minus_array[0][ei]);
    }
    
    return beta_minus;
}


/*
 * Perform WENO reconstruction.
 */
void
ConvectiveFluxReconstructorWENO_CU6_M2_LLF::performWENOReconstruction(
   std::vector<double>& G_plus,
   std::vector<double>& G_minus,
   const boost::multi_array<double, 2>& G_plus_array,
   const boost::multi_array<double, 2>& G_minus_array,
   const DIRECTION& direction)
{
    const double epsilon = d_constant_epsilon;
    const double Chi = d_constant_Chi;
    
    const double* const grid_spacing = d_grid_geometry->getDx();
    double dx = 0.0;
    if (direction == X_DIRECTION)
    {
        dx = grid_spacing[0];
    }
    else if (direction == Y_DIRECTION)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(d_dim.getValue() >= 2.0);
#endif
        dx = grid_spacing[1];
    }
    else if (direction == Z_DIRECTION)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(d_dim.getValue() >= 3.0);
#endif
        dx = grid_spacing[2];
    }
    
    /*
     * Compute the beta_plus's and beta_minus's.
     */
    
    boost::multi_array<double, 2> beta_plus = computeBetaPlus(G_plus_array);
    boost::multi_array<double, 2> beta_minus = computeBetaMinus(G_minus_array);
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        /*
         * Do WENO reconstruction on current characteristic variables.
         */
        
        std::vector<double> G_plus_r;
        std::vector<double> G_minus_r;
        
        // Do the reconstruction from different stencils.
        for (int r = 0; r < 4; r++)
        {
            G_plus_r.push_back(0.0);
            G_minus_r.push_back(0.0);
            
            for (int m = r; m < 3 + r; m++)
            {
                G_plus_r[r] += d_weights_c[r][m - r]*G_plus_array[m][ei];
                G_minus_r[r] += d_weights_c[r][m - r]*G_minus_array[6 - m - 1][ei];
            }
        }
        
        /*
         * Compute G_plus of the current characteristic variable.
         */
        
        // Compute the reference smoothness indicators tau_6_plus.
        const double beta_plus_average = 1.0/6*(beta_plus[ei][0] + beta_plus[ei][2] +
            4*beta_plus[ei][1]);
        
        const double tau_6_plus = beta_plus[ei][3] - beta_plus_average;
        
        // Compute the weights alpha_plus.
        std::vector<double> alpha_plus;
        for (int r = 0; r < 4; r++)
        {
            alpha_plus.push_back(
                d_weights_d[r]*pow(d_constant_C + tau_6_plus/(beta_plus[ei][r] + epsilon*dx*dx)*
                    (beta_plus_average + Chi*dx*dx)/(beta_plus[ei][r] + Chi*dx*dx), d_constant_q));
        }
        
        // Sum up the weights alpha_plus.
        double alpha_plus_sum = 0.0;
        for (int r = 0; r < 4; r++)
        {
            alpha_plus_sum += alpha_plus[r];
        }
        
        // Compute the nonlinear weights omega_plus.
        std::vector<double> omega_plus;
        for (int r = 0; r < 4; r++)
        {
            omega_plus.push_back(alpha_plus[r]/alpha_plus_sum);
        }
        
        // Compute the G_plus.
        G_plus.push_back(0.0);
        for (int r = 0; r < 4; r++)
        {
            G_plus[ei] += omega_plus[r]*G_plus_r[r];
        }
        
        /*
         * Compute G_minus of the current characteristic variable.
         */
        
        // Compute the reference smoothness indicators tau_6_plus.
        const double beta_minus_average = 1.0/6*(beta_minus[ei][0] + beta_minus[ei][2] +
            4*beta_minus[ei][1]);
        
        const double tau_6_minus = beta_minus[ei][3] - beta_minus_average;
        
        // Compute the weights alpha_minus.
        std::vector<double> alpha_minus;
        for (int r = 0; r < 4; r++)
        {
            alpha_minus.push_back(
                d_weights_d[r]*pow(d_constant_C + tau_6_minus/(beta_minus[ei][r] + epsilon*dx*dx)*
                    (beta_minus_average + Chi*dx*dx)/(beta_minus[ei][r] + Chi*dx*dx), d_constant_q));
        }
        
        // Sum up the weights alpha_minus.
        double alpha_minus_sum = 0.0;
        for (int r = 0; r < 4; r++)
        {
            alpha_minus_sum += alpha_minus[r];
        }
        
        // Compute the nonlinear weights omega_minus.
        std::vector<double> omega_minus;
        for (int r = 0; r < 4; r++)
        {
            omega_minus.push_back(alpha_minus[r]/alpha_minus_sum);
        }
        
        // Compute the G_minus.
        G_minus.push_back(0.0);
        for (int r = 0; r < 4; r++)
        {
            G_minus[ei] += omega_minus[r]*G_minus_r[r];
        }
    }
}
