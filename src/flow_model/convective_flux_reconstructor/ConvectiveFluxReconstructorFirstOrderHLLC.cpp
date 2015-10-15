#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorFirstOrderHLLC.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

ConvectiveFluxReconstructorFirstOrderHLLC::ConvectiveFluxReconstructorFirstOrderHLLC(
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
            shock_capturing_scheme_db),
        d_Riemann_solver(
            "HLLC Riemann solver",
            d_dim,
            d_num_eqn,
            d_num_species,
            d_equation_of_state)
{
    d_num_conv_ghosts = hier::IntVector::getOne(d_dim);
}

/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorFirstOrderHLLC::printClassData(
    std::ostream& os) const
{
    os << "\nConvectiveFluxReconstructorFirstOrderHLLC::printClassData..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorFirstOrderHLLC: this = "
       << (ConvectiveFluxReconstructorFirstOrderHLLC *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_set_variables = "
       << d_set_variables
       << std::endl;
    
    os << std::endl;
    
    os << "End of ConvectiveFluxReconstructorFirstOrderHLLC::printClassData"
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorFirstOrderHLLC::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putString("d_shock_capturing_scheme", "FIRST_ORDER_HLLC");
}


/*
 * Compute the convective fluxes and sources due to hyperbolization
 * of the equations.
 */
void
ConvectiveFluxReconstructorFirstOrderHLLC::computeConvectiveFluxesAndSources(
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
                
                if (d_dim == tbox::Dimension(1))
                {
                    // Get the arrays of time-dependent variables.
                    double* rho   = density->getPointer(0);
                    double* rho_u = momentum->getPointer(0);
                    double* E     = total_energy->getPointer(0);
                    
                    // Compute the fluxes in the x direction.
                    for (int i = 0; i < interior_dims[0] + 1; i++)
                    {
                        // Compute the indices of left cell, right cell and face.
                        const int idx_cell_L = i - 1 + d_num_ghosts[0];
                        
                        const int idx_cell_R = i + d_num_ghosts[0];
                        
                        const int idx_face_x = i;
                        
                        std::vector<double*> F_x_ptr;
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            F_x_ptr.push_back(&(convective_flux->getPointer(0, ei)[idx_face_x]));
                        }
                        
                        std::vector<const double*> m_L_ptr;
                        m_L_ptr.push_back(&rho_u[idx_cell_L]);
                        
                        std::vector<const double*> m_R_ptr;
                        m_R_ptr.push_back(&rho_u[idx_cell_R]);
                        
                        d_Riemann_solver.computeIntercellFluxForSingleSpecies(
                            F_x_ptr,
                            &rho[idx_cell_L],
                            &rho[idx_cell_R],
                            m_L_ptr,
                            m_R_ptr,
                            &E[idx_cell_L],
                            &E[idx_cell_R],
                            X_DIRECTION);
                        
                        // Mulitply fluxes by dt.
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            *F_x_ptr[ei] *= dt;
                        }
                    }
                }  // if (d_dim == tbox::Dimension(1))
                else if (d_dim == tbox::Dimension(2))
                {
                    // Get the arrays of time-dependent variables.
                    double* rho   = density->getPointer(0);
                    double* rho_u = momentum->getPointer(0);
                    double* rho_v = momentum->getPointer(1);
                    double* E     = total_energy->getPointer(0);
                    
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
                            
                            std::vector<double*> F_x_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_x_ptr.push_back(&(convective_flux->getPointer(0, ei)[idx_face_x]));
                            }
                            
                            std::vector<const double*> m_L_ptr;
                            m_L_ptr.push_back(&rho_u[idx_cell_L]);
                            m_L_ptr.push_back(&rho_v[idx_cell_L]);
                            
                            std::vector<const double*> m_R_ptr;
                            m_R_ptr.push_back(&rho_u[idx_cell_R]);
                            m_R_ptr.push_back(&rho_v[idx_cell_R]);
                            
                            d_Riemann_solver.computeIntercellFluxForSingleSpecies(
                                F_x_ptr,
                                &rho[idx_cell_L],
                                &rho[idx_cell_R],
                                m_L_ptr,
                                m_R_ptr,
                                &E[idx_cell_L],
                                &E[idx_cell_R],
                                X_DIRECTION);
                            
                            // Mulitply fluxes by dt.
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                *F_x_ptr[ei] *= dt;
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
                            
                            std::vector<double*> F_y_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_y_ptr.push_back(&(convective_flux->getPointer(1, ei)[idx_face_y]));
                            }
                            
                            std::vector<const double*> m_B_ptr;
                            m_B_ptr.push_back(&rho_u[idx_cell_B]);
                            m_B_ptr.push_back(&rho_v[idx_cell_B]);
                            
                            std::vector<const double*> m_T_ptr;
                            m_T_ptr.push_back(&rho_u[idx_cell_T]);
                            m_T_ptr.push_back(&rho_v[idx_cell_T]);
                            
                            d_Riemann_solver.computeIntercellFluxForSingleSpecies(
                                F_y_ptr,
                                &rho[idx_cell_B],
                                &rho[idx_cell_T],
                                m_B_ptr,
                                m_T_ptr,
                                &E[idx_cell_B],
                                &E[idx_cell_T],
                                Y_DIRECTION);
                            
                            // Mulitply fluxes by dt.
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                *F_y_ptr[ei] *= dt;
                            }
                        }
                    }
                } // if (d_dim == tbox::Dimension(2))
                else if (d_dim == tbox::Dimension(3))
                {
                    // Get the arrays of time-dependent variables.
                    double* rho   = density->getPointer(0);
                    double* rho_u = momentum->getPointer(0);
                    double* rho_v = momentum->getPointer(1);
                    double* rho_w = momentum->getPointer(2);
                    double* E     = total_energy->getPointer(0);
                    
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
                                
                                std::vector<double*> F_x_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_x_ptr.push_back(&(convective_flux->getPointer(0, ei)[idx_face_x]));
                                }
                                
                                std::vector<const double*> m_L_ptr;
                                m_L_ptr.push_back(&rho_u[idx_cell_L]);
                                m_L_ptr.push_back(&rho_v[idx_cell_L]);
                                m_L_ptr.push_back(&rho_w[idx_cell_L]);
                                
                                std::vector<const double*> m_R_ptr;
                                m_R_ptr.push_back(&rho_u[idx_cell_R]);
                                m_R_ptr.push_back(&rho_v[idx_cell_R]);
                                m_R_ptr.push_back(&rho_w[idx_cell_R]);
                                
                                d_Riemann_solver.computeIntercellFluxForSingleSpecies(
                                    F_x_ptr,
                                    &rho[idx_cell_L],
                                    &rho[idx_cell_R],
                                    m_L_ptr,
                                    m_R_ptr,
                                    &E[idx_cell_L],
                                    &E[idx_cell_R],
                                    X_DIRECTION);
                                
                                // Mulitply fluxes by dt.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    *F_x_ptr[ei] *= dt;
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
                                
                                std::vector<double*> F_y_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_y_ptr.push_back(&(convective_flux->getPointer(1, ei)[idx_face_y]));
                                }
                                
                                std::vector<const double*> m_B_ptr;
                                m_B_ptr.push_back(&rho_u[idx_cell_B]);
                                m_B_ptr.push_back(&rho_v[idx_cell_B]);
                                m_B_ptr.push_back(&rho_w[idx_cell_B]);
                                
                                std::vector<const double*> m_T_ptr;
                                m_T_ptr.push_back(&rho_u[idx_cell_T]);
                                m_T_ptr.push_back(&rho_v[idx_cell_T]);
                                m_T_ptr.push_back(&rho_w[idx_cell_T]);
                                
                                d_Riemann_solver.computeIntercellFluxForSingleSpecies(
                                    F_y_ptr,
                                    &rho[idx_cell_B],
                                    &rho[idx_cell_T],
                                    m_B_ptr,
                                    m_T_ptr,
                                    &E[idx_cell_B],
                                    &E[idx_cell_T],
                                    Y_DIRECTION);
                                
                                // Mulitply fluxes by dt.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    *F_y_ptr[ei] *= dt;
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
                                
                                std::vector<double*> F_z_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_z_ptr.push_back(&(convective_flux->getPointer(2, ei)[idx_face_z]));
                                }
                                
                                std::vector<const double*> m_B_ptr;
                                m_B_ptr.push_back(&rho_u[idx_cell_B]);
                                m_B_ptr.push_back(&rho_v[idx_cell_B]);
                                m_B_ptr.push_back(&rho_w[idx_cell_B]);
                                
                                std::vector<const double*> m_F_ptr;
                                m_F_ptr.push_back(&rho_u[idx_cell_F]);
                                m_F_ptr.push_back(&rho_v[idx_cell_F]);
                                m_F_ptr.push_back(&rho_w[idx_cell_F]);
                                
                                d_Riemann_solver.computeIntercellFluxForSingleSpecies(
                                    F_z_ptr,
                                    &rho[idx_cell_B],
                                    &rho[idx_cell_F],
                                    m_B_ptr,
                                    m_F_ptr,
                                    &E[idx_cell_B],
                                    &E[idx_cell_F],
                                    Z_DIRECTION);
                                
                                // Mulitply fluxes by dt.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    *F_z_ptr[ei] *= dt;
                                }
                            }
                        }
                    }
                } // if (d_dim == tbox::Dimension(3))
                
                break;
            }
            case FOUR_EQN_SHYUE:
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
                
                boost::shared_ptr<pdat::CellData<double> > mass_fraction(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_mass_fraction, data_context)));
                
                boost::shared_ptr<pdat::FaceData<double> > convective_flux(
                    BOOST_CAST<pdat::FaceData<double>, hier::PatchData>(
                        patch.getPatchData(d_convective_flux, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > source(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_source, data_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(density);
                TBOX_ASSERT(momentum);
                TBOX_ASSERT(total_energy);
                TBOX_ASSERT(mass_fraction);
                TBOX_ASSERT(convective_flux);
                TBOX_ASSERT(source);
                
                TBOX_ASSERT(density->getGhostCellWidth() == d_num_ghosts);
                TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
                TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
                TBOX_ASSERT(mass_fraction->getGhostCellWidth() == d_num_ghosts);
                TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
                TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
                
                // Allocate temporary patch data.
                boost::shared_ptr<pdat::FaceData<double> > velocity_intercell(
                    new pdat::FaceData<double>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
                
                if (d_dim == tbox::Dimension(1))
                {
                    // Get the arrays of time-dependent variables.
                    double* rho   = density->getPointer(0);
                    double* rho_u = momentum->getPointer(0);
                    double* E     = total_energy->getPointer(0);
                    std::vector<double*> Y;
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Y.push_back(mass_fraction->getPointer(si));
                    }
                    
                    // Compute the fluxes in the x direction and velocity components at the face
                    // normal to the x direction.
                    for (int i = 0; i < interior_dims[0] + 1; i++)
                    {
                        // Compute the indices of left cell, right cell and face.
                        const int idx_cell_L = i - 1 + d_num_ghosts[0];
                        
                        const int idx_cell_R = i + d_num_ghosts[0];
                        
                        const int idx_face_x = i;
                        
                        std::vector<double*> F_x_ptr;
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            F_x_ptr.push_back(&(convective_flux->getPointer(0, ei)[idx_face_x]));
                        }
                        
                        std::vector<double*> vel_x_ptr;
                        vel_x_ptr.push_back(&(velocity_intercell->getPointer(0, 0)[idx_face_x]));
                        
                        std::vector<const double*> m_L_ptr;
                        m_L_ptr.push_back(&rho_u[idx_cell_L]);
                        
                        std::vector<const double*> m_R_ptr;
                        m_R_ptr.push_back(&rho_u[idx_cell_R]);
                        
                        std::vector<const double*> Y_L_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Y_L_ptr.push_back(
                                &Y[si][idx_cell_L]);
                        }
                        
                        std::vector<const double*> Y_R_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Y_R_ptr.push_back(
                                &Y[si][idx_cell_R]);
                        }
                        
                        d_Riemann_solver.computeIntercellFluxAndVelocityForFourEqnShyue(
                            F_x_ptr,
                            vel_x_ptr,
                            &rho[idx_cell_L],
                            &rho[idx_cell_R],
                            m_L_ptr,
                            m_R_ptr,
                            &E[idx_cell_L],
                            &E[idx_cell_R],
                            Y_L_ptr,
                            Y_R_ptr,
                            X_DIRECTION);
                        
                        // Mulitply fluxes by dt.
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            *F_x_ptr[ei] *= dt;
                        }
                    }
                    
                    // Compute the source.
                    // Get the grid spacing.
                    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                            patch.getPatchGeometry()));
                    
                    const double* dx = patch_geom->getDx();
                    
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        double *S = source->getPointer(d_num_eqn - (d_num_species - 1) + si);
                        
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute the indices of cell and faces. 
                            const int idx_cell_wghost = i + d_num_ghosts[0];
                            
                            const int idx_cell_nghost = i;
                            
                            const int idx_face_x_L = i;
                            
                            const int idx_face_x_R = i + 1;
                            
                            const int& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                            const int& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                            
                            S[idx_cell_nghost] += dt*Y[si][idx_cell_wghost]*(u_R - u_L)/dx[0];
                        }
                    }
                } // if (d_dim == tbox::Dimension(1))
                else if (d_dim == tbox::Dimension(2))
                {
                    // Get the arrays of time-dependent variables.
                    double* rho   = density->getPointer(0);
                    double* rho_u = momentum->getPointer(0);
                    double* rho_v = momentum->getPointer(1);
                    double* E     = total_energy->getPointer(0);
                    std::vector<double*> Y;
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Y.push_back(mass_fraction->getPointer(si));
                    }
                    
                    // Compute the fluxes in the x direction and velocity components at the face
                    // normal to the x direction.
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
                            
                            std::vector<double*> F_x_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_x_ptr.push_back(&(convective_flux->getPointer(0, ei)[idx_face_x]));
                            }
                            
                            std::vector<double*> vel_x_ptr;
                            for (int vi = 0; vi < 2; vi++)
                            {
                                vel_x_ptr.push_back(&(velocity_intercell->getPointer(0, vi)[idx_face_x]));
                            }
                            
                            std::vector<const double*> m_L_ptr;
                            m_L_ptr.push_back(&rho_u[idx_cell_L]);
                            m_L_ptr.push_back(&rho_v[idx_cell_L]);
                            
                            std::vector<const double*> m_R_ptr;
                            m_R_ptr.push_back(&rho_u[idx_cell_R]);
                            m_R_ptr.push_back(&rho_v[idx_cell_R]);
                            
                            std::vector<const double*> Y_L_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Y_L_ptr.push_back(
                                    &Y[si][idx_cell_L]);
                            }
                            
                            std::vector<const double*> Y_R_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Y_R_ptr.push_back(
                                    &Y[si][idx_cell_R]);
                            }
                            
                            d_Riemann_solver.computeIntercellFluxAndVelocityForFourEqnShyue(
                                F_x_ptr,
                                vel_x_ptr,
                                &rho[idx_cell_L],
                                &rho[idx_cell_R],
                                m_L_ptr,
                                m_R_ptr,
                                &E[idx_cell_L],
                                &E[idx_cell_R],
                                Y_L_ptr,
                                Y_R_ptr,
                                X_DIRECTION);
                            
                            // Mulitply fluxes by dt.
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                *F_x_ptr[ei] *= dt;
                            }
                        }
                    }
                    
                    // Compute the fluxes in the y direction and velocity components at the face
                    // normal to the y direction.
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
                            
                            std::vector<double*> F_y_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_y_ptr.push_back(&(convective_flux->getPointer(1, ei)[idx_face_y]));
                            }
                            
                            std::vector<double*> vel_y_ptr;
                            for (int vi = 0; vi < 2; vi++)
                            {
                                vel_y_ptr.push_back(&(velocity_intercell->getPointer(1, vi)[idx_face_y]));
                            }
                            
                            std::vector<const double*> m_B_ptr;
                            m_B_ptr.push_back(&rho_u[idx_cell_B]);
                            m_B_ptr.push_back(&rho_v[idx_cell_B]);
                            
                            std::vector<const double*> m_T_ptr;
                            m_T_ptr.push_back(&rho_u[idx_cell_T]);
                            m_T_ptr.push_back(&rho_v[idx_cell_T]);
                            
                            std::vector<const double*> Y_B_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Y_B_ptr.push_back(
                                    &Y[si][idx_cell_B]);
                            }
                            
                            std::vector<const double*> Y_T_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Y_T_ptr.push_back(
                                    &Y[si][idx_cell_T]);
                            }
                            
                            d_Riemann_solver.computeIntercellFluxAndVelocityForFourEqnShyue(
                                F_y_ptr,
                                vel_y_ptr,
                                &rho[idx_cell_B],
                                &rho[idx_cell_T],
                                m_B_ptr,
                                m_T_ptr,
                                &E[idx_cell_B],
                                &E[idx_cell_T],
                                Y_B_ptr,
                                Y_T_ptr,
                                Y_DIRECTION);
                            
                            // Mulitply fluxes by dt.
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                *F_y_ptr[ei] *= dt;
                            }
                        }
                    }
                    
                    // Compute the source.
                    // Get the grid spacing.
                    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                            patch.getPatchGeometry()));
                    
                    const double* dx = patch_geom->getDx();
                    
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        double *S = source->getPointer(d_num_eqn - (d_num_species - 1) + si);
                        
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0]; i++)
                            {
                                // Compute the indices of cell and faces. 
                                const int idx_cell_wghost = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_nghost = i + j*interior_dims[0];
                                
                                const int idx_face_x_L = i +
                                    j*(interior_dims[0] + 1);
                                
                                const int idx_face_x_R = (i + 1) +
                                    j*(interior_dims[0] + 1);
                                
                                const int idx_face_y_B = j +
                                    i*(interior_dims[1] + 1);
                                
                                const int idx_face_y_T = (j + 1) +
                                    i*(interior_dims[1] + 1);
                                
                                const int& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                                const int& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                                
                                const int& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                                const int& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                                
                                S[idx_cell_nghost] += dt*Y[si][idx_cell_wghost]*((u_R - u_L)/dx[0] + (v_T - v_B)/dx[1]);
                            }
                        }
                    }
                } // if (d_dim == tbox::Dimension(2))
                else if (d_dim == tbox::Dimension(3))
                {
                    // Get the arrays of time-dependent variables.
                    double* rho   = density->getPointer(0);
                    double* rho_u = momentum->getPointer(0);
                    double* rho_v = momentum->getPointer(1);
                    double* rho_w = momentum->getPointer(2);
                    double* E     = total_energy->getPointer(0);
                    std::vector<double*> Y;
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Y.push_back(mass_fraction->getPointer(si));
                    }
                    
                    // Compute the fluxes in the x direction and velocity components at the face
                    // normal to the x direction.
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
                                
                                std::vector<double*> F_x_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_x_ptr.push_back(&(convective_flux->getPointer(0, ei)[idx_face_x]));
                                }
                                
                                std::vector<double*> vel_x_ptr;
                                for (int vi = 0; vi < 3; vi++)
                                {
                                    vel_x_ptr.push_back(&(velocity_intercell->getPointer(0, vi)[idx_face_x]));
                                }
                                
                                std::vector<const double*> m_L_ptr;
                                m_L_ptr.push_back(&rho_u[idx_cell_L]);
                                m_L_ptr.push_back(&rho_v[idx_cell_L]);
                                m_L_ptr.push_back(&rho_w[idx_cell_L]);
                                
                                std::vector<const double*> m_R_ptr;
                                m_R_ptr.push_back(&rho_u[idx_cell_R]);
                                m_R_ptr.push_back(&rho_v[idx_cell_R]);
                                m_R_ptr.push_back(&rho_w[idx_cell_R]);
                                
                                std::vector<const double*> Y_L_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_L_ptr.push_back(
                                        &Y[si][idx_cell_L]);
                                }
                                
                                std::vector<const double*> Y_R_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_R_ptr.push_back(
                                        &Y[si][idx_cell_R]);
                                }
                                
                                d_Riemann_solver.computeIntercellFluxAndVelocityForFourEqnShyue(
                                    F_x_ptr,
                                    vel_x_ptr,
                                    &rho[idx_cell_L],
                                    &rho[idx_cell_R],
                                    m_L_ptr,
                                    m_R_ptr,
                                    &E[idx_cell_L],
                                    &E[idx_cell_R],
                                    Y_L_ptr,
                                    Y_R_ptr,
                                    X_DIRECTION);
                                
                                // Mulitply fluxes by dt.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    *F_x_ptr[ei] *= dt;
                                }
                            }
                        }
                    }
                    
                    // Compute the fluxes in the y direction and velocity components at the face
                    // normal to the y direction.
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
                                
                                std::vector<double*> F_y_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_y_ptr.push_back(&(convective_flux->getPointer(1, ei)[idx_face_y]));
                                }
                                
                                std::vector<double*> vel_y_ptr;
                                for (int vi = 0; vi < 3; vi++)
                                {
                                    vel_y_ptr.push_back(&(velocity_intercell->getPointer(1, vi)[idx_face_y]));
                                }
                                
                                std::vector<const double*> m_B_ptr;
                                m_B_ptr.push_back(&rho_u[idx_cell_B]);
                                m_B_ptr.push_back(&rho_v[idx_cell_B]);
                                m_B_ptr.push_back(&rho_w[idx_cell_B]);
                                
                                std::vector<const double*> m_T_ptr;
                                m_T_ptr.push_back(&rho_u[idx_cell_T]);
                                m_T_ptr.push_back(&rho_v[idx_cell_T]);
                                m_T_ptr.push_back(&rho_w[idx_cell_T]);
                                
                                std::vector<const double*> Y_B_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_B_ptr.push_back(
                                        &Y[si][idx_cell_B]);
                                }
                                
                                std::vector<const double*> Y_T_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_T_ptr.push_back(
                                        &Y[si][idx_cell_T]);
                                }
                                
                                d_Riemann_solver.computeIntercellFluxAndVelocityForFourEqnShyue(
                                    F_y_ptr,
                                    vel_y_ptr,
                                    &rho[idx_cell_B],
                                    &rho[idx_cell_T],
                                    m_B_ptr,
                                    m_T_ptr,
                                    &E[idx_cell_B],
                                    &E[idx_cell_T],
                                    Y_B_ptr,
                                    Y_T_ptr,
                                    Y_DIRECTION);
                                
                                // Mulitply fluxes by dt.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    *F_y_ptr[ei] *= dt;
                                }
                            }
                        }
                    }
                    
                    // Compute the fluxes in the z direction and velocity components at the face
                    // normal to the z direction.
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
                                
                                std::vector<double*> F_z_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_z_ptr.push_back(&(convective_flux->getPointer(2, ei)[idx_face_z]));
                                }
                                
                                std::vector<double*> vel_z_ptr;
                                for (int vi = 0; vi < 3; vi++)
                                {
                                    vel_z_ptr.push_back(&(velocity_intercell->getPointer(2, vi)[idx_face_z]));
                                }
                                
                                std::vector<const double*> m_B_ptr;
                                m_B_ptr.push_back(&rho_u[idx_cell_B]);
                                m_B_ptr.push_back(&rho_v[idx_cell_B]);
                                m_B_ptr.push_back(&rho_w[idx_cell_B]);
                                
                                std::vector<const double*> m_F_ptr;
                                m_F_ptr.push_back(&rho_u[idx_cell_F]);
                                m_F_ptr.push_back(&rho_v[idx_cell_F]);
                                m_F_ptr.push_back(&rho_w[idx_cell_F]);
                                
                                std::vector<const double*> Y_B_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_B_ptr.push_back(
                                        &Y[si][idx_cell_B]);
                                }
                                
                                std::vector<const double*> Y_F_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_F_ptr.push_back(
                                        &Y[si][idx_cell_F]);
                                }
                                
                                d_Riemann_solver.computeIntercellFluxAndVelocityForFourEqnShyue(
                                    F_z_ptr,
                                    vel_z_ptr,
                                    &rho[idx_cell_B],
                                    &rho[idx_cell_F],
                                    m_B_ptr,
                                    m_F_ptr,
                                    &E[idx_cell_B],
                                    &E[idx_cell_F],
                                    Y_B_ptr,
                                    Y_F_ptr,
                                    Z_DIRECTION);
                                
                                // Mulitply fluxes by dt.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    *F_z_ptr[ei] *= dt;
                                }
                            }
                        }
                    }
                    
                    // Compute the source.
                    // Get the grid spacing.
                    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                            patch.getPatchGeometry()));
                    
                    const double* dx = patch_geom->getDx();
                    
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        double *S = source->getPointer(d_num_eqn - (d_num_species - 1) + si);
                        
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = 0; j < interior_dims[1]; j++)
                            {
                                for (int i = 0; i < interior_dims[0]; i++)
                                {
                                    // Compute the indices of cell and faces. 
                                    const int idx_cell_wghost = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_nghost = i +
                                        j*interior_dims[0] +
                                        k*interior_dims[0]*interior_dims[1];
                                    
                                    const int idx_face_x_L = i +
                                        j*(interior_dims[0] + 1) +
                                        k*(interior_dims[0] + 1)*interior_dims[1];
                                    
                                    const int idx_face_x_R = (i + 1) +
                                        j*(interior_dims[0] + 1) +
                                        k*(interior_dims[0] + 1)*interior_dims[1];
                                    
                                    const int idx_face_y_B = j +
                                        k*(interior_dims[1] + 1) +
                                        i*(interior_dims[1] + 1)*interior_dims[2];
                                    
                                    const int idx_face_y_T = (j + 1) +
                                        k*(interior_dims[1] + 1) +
                                        i*(interior_dims[1] + 1)*interior_dims[2];
                                    
                                    const int idx_face_z_B = k +
                                        i*(interior_dims[2] + 1) +
                                        j*(interior_dims[2] + 1)*interior_dims[0];
                                    
                                    const int idx_face_z_F = (k + 1) +
                                        i*(interior_dims[2] + 1) +
                                        j*(interior_dims[2] + 1)*interior_dims[0];
                                    
                                    const int& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                                    const int& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                                    
                                    const int& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                                    const int& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                                    
                                    const int& w_B = velocity_intercell->getPointer(2, 2)[idx_face_z_B];
                                    const int& w_F = velocity_intercell->getPointer(2, 2)[idx_face_z_F];
                                    
                                    S[idx_cell_nghost] += dt*Y[si][idx_cell_wghost]*((u_R - u_L)/dx[0] +
                                        (v_T - v_B)/dx[1] + (w_F - w_B)/dx[2]);
                                }
                            }
                        }
                    }
                } // if (d_dim == tbox::Dimension(3))
                
                break;
            }
            case FIVE_EQN_ALLAIRE:
            {
                // Get the cell data of the time-dependent variables.
                boost::shared_ptr<pdat::CellData<double> > partial_density(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_partial_density, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > momentum(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_momentum, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > total_energy(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_total_energy, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > volume_fraction(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_volume_fraction, data_context)));
                
                boost::shared_ptr<pdat::FaceData<double> > convective_flux(
                    BOOST_CAST<pdat::FaceData<double>, hier::PatchData>(
                        patch.getPatchData(d_convective_flux, data_context)));
                
                boost::shared_ptr<pdat::CellData<double> > source(
                    BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch.getPatchData(d_source, data_context)));
                
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(partial_density);
        TBOX_ASSERT(momentum);
        TBOX_ASSERT(total_energy);
        TBOX_ASSERT(volume_fraction);
        TBOX_ASSERT(convective_flux);
        TBOX_ASSERT(source);
        
        TBOX_ASSERT(partial_density->getGhostCellWidth() == d_num_ghosts);
        TBOX_ASSERT(momentum->getGhostCellWidth() == d_num_ghosts);
        TBOX_ASSERT(total_energy->getGhostCellWidth() == d_num_ghosts);
        TBOX_ASSERT(volume_fraction->getGhostCellWidth() == d_num_ghosts);
        TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
        TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
                
                // Allocate temporary patch data.
                boost::shared_ptr<pdat::FaceData<double> > velocity_intercell(
                    new pdat::FaceData<double>(interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
                
                if (d_dim == tbox::Dimension(1))
                {
                    // Get the arrays of time-dependent variables.
                    std::vector<double*> Z_rho;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Z_rho.push_back(partial_density->getPointer(si));
                    }  
                    double* rho_u = momentum->getPointer(0);
                    double* E     = total_energy->getPointer(0);
                    std::vector<double*> Z;
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Z.push_back(volume_fraction->getPointer(si));
                    }
                    
                    // Compute the fluxes in the x direction and velocity components at the face
                    // normal to the x direction.
                    for (int i = 0; i < interior_dims[0] + 1; i++)
                    {
                        // Compute the indices of left cell, right cell and face.
                        const int idx_cell_L = i - 1 + d_num_ghosts[0];
                        
                        const int idx_cell_R = i + d_num_ghosts[0];
                        
                        const int idx_face_x = i;
                        
                        std::vector<double*> F_x_ptr;
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            F_x_ptr.push_back(&(convective_flux->getPointer(0, ei)[idx_face_x]));
                        }
                        
                        std::vector<double*> vel_x_ptr;
                        vel_x_ptr.push_back(&(velocity_intercell->getPointer(0, 0)[idx_face_x]));
                        
                        std::vector<const double*> Z_rho_L_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho_L_ptr.push_back(
                                &Z_rho[si][idx_cell_L]);
                        }
                        
                        std::vector<const double*> Z_rho_R_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho_R_ptr.push_back(
                                &Z_rho[si][idx_cell_R]);
                        }
                        
                        std::vector<const double*> m_L_ptr;
                        m_L_ptr.push_back(&rho_u[idx_cell_L]);
                        
                        std::vector<const double*> m_R_ptr;
                        m_R_ptr.push_back(&rho_u[idx_cell_R]);
                        
                        std::vector<const double*> Z_L_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Z_L_ptr.push_back(
                                &Z[si][idx_cell_L]);
                        }
                        
                        std::vector<const double*> Z_R_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Z_R_ptr.push_back(
                                &Z[si][idx_cell_R]);
                        }
                        
                        d_Riemann_solver.computeIntercellFluxAndVelocityForFiveEqnAllaire(
                            F_x_ptr,
                            vel_x_ptr,
                            Z_rho_L_ptr,
                            Z_rho_R_ptr,
                            m_L_ptr,
                            m_R_ptr,
                            &E[idx_cell_L],
                            &E[idx_cell_R],
                            Z_L_ptr,
                            Z_R_ptr,
                            X_DIRECTION);
                        
                        // Mulitply fluxes by dt.
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            *F_x_ptr[ei] *= dt;
                        }
                    }
                    
                    // Compute the source.
                    // Get the grid spacing.
                    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                            patch.getPatchGeometry()));
                    
                    const double* dx = patch_geom->getDx();
                    
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        double *S = source->getPointer(d_num_eqn - (d_num_species - 1) + si);
                        
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute the indices of cell and faces. 
                            const int idx_cell_wghost = i + d_num_ghosts[0];
                            
                            const int idx_cell_nghost = i;
                            
                            const int idx_face_x_L = i;
                            
                            const int idx_face_x_R = i + 1;
                            
                            const int& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                            const int& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                            
                            S[idx_cell_nghost] += dt*Z[si][idx_cell_wghost]*(u_R - u_L)/dx[0];
                        }
                    }
                }  // if (d_dim == tbox::Dimension(1))
                else if (d_dim == tbox::Dimension(2))
                {
                    // Get the arrays of time-dependent variables.
                    std::vector<double*> Z_rho;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Z_rho.push_back(partial_density->getPointer(si));
                    }  
                    double* rho_u = momentum->getPointer(0);
                    double* rho_v = momentum->getPointer(1);
                    double* E     = total_energy->getPointer(0);
                    std::vector<double*> Z;
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Z.push_back(volume_fraction->getPointer(si));
                    }
                    
                    // Compute the fluxes in the x direction and velocity components at the face
                    // normal to the x direction.
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
                            
                            std::vector<double*> F_x_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_x_ptr.push_back(&(convective_flux->getPointer(0, ei)[idx_face_x]));
                            }
                            
                            std::vector<double*> vel_x_ptr;
                            for (int vi = 0; vi < 2; vi++)
                            {
                                vel_x_ptr.push_back(&(velocity_intercell->getPointer(0, vi)[idx_face_x]));
                            }
                            
                            std::vector<const double*> Z_rho_L_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_L_ptr.push_back(
                                    &Z_rho[si][idx_cell_L]);
                            }
                            
                            std::vector<const double*> Z_rho_R_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_R_ptr.push_back(
                                    &Z_rho[si][idx_cell_R]);
                            }
                            
                            std::vector<const double*> m_L_ptr;
                            m_L_ptr.push_back(&rho_u[idx_cell_L]);
                            m_L_ptr.push_back(&rho_v[idx_cell_L]);
                            
                            std::vector<const double*> m_R_ptr;
                            m_R_ptr.push_back(&rho_u[idx_cell_R]);
                            m_R_ptr.push_back(&rho_v[idx_cell_R]);
                            
                            std::vector<const double*> Z_L_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Z_L_ptr.push_back(
                                    &Z[si][idx_cell_L]);
                            }
                            
                            std::vector<const double*> Z_R_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Z_R_ptr.push_back(
                                    &Z[si][idx_cell_R]);
                            }
                            
                            d_Riemann_solver.computeIntercellFluxAndVelocityForFiveEqnAllaire(
                                F_x_ptr,
                                vel_x_ptr,
                                Z_rho_L_ptr,
                                Z_rho_R_ptr,
                                m_L_ptr,
                                m_R_ptr,
                                &E[idx_cell_L],
                                &E[idx_cell_R],
                                Z_L_ptr,
                                Z_R_ptr,
                                X_DIRECTION);
                            
                            // Mulitply fluxes by dt.
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                *F_x_ptr[ei] *= dt;
                            }
                        }
                    }
                    
                    // Compute the fluxes in the y direction and velocity components at the face
                    // normal to the y direction.
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
                            
                            std::vector<double*> F_y_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_y_ptr.push_back(&(convective_flux->getPointer(1, ei)[idx_face_y]));
                            }
                            
                            std::vector<double*> vel_y_ptr;
                            for (int vi = 0; vi < 2; vi++)
                            {
                                vel_y_ptr.push_back(&(velocity_intercell->getPointer(1, vi)[idx_face_y]));
                            }
                            
                            std::vector<const double*> Z_rho_B_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_B_ptr.push_back(
                                    &Z_rho[si][idx_cell_B]);
                            }
                            
                            std::vector<const double*> Z_rho_T_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_T_ptr.push_back(
                                    &Z_rho[si][idx_cell_T]);
                            }
                            
                            std::vector<const double*> m_B_ptr;
                            m_B_ptr.push_back(&rho_u[idx_cell_B]);
                            m_B_ptr.push_back(&rho_v[idx_cell_B]);
                            
                            std::vector<const double*> m_T_ptr;
                            m_T_ptr.push_back(&rho_u[idx_cell_T]);
                            m_T_ptr.push_back(&rho_v[idx_cell_T]);
                            
                            std::vector<const double*> Z_B_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Z_B_ptr.push_back(
                                    &Z[si][idx_cell_B]);
                            }
                            
                            std::vector<const double*> Z_T_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Z_T_ptr.push_back(
                                    &Z[si][idx_cell_T]);
                            }
                            
                            d_Riemann_solver.computeIntercellFluxAndVelocityForFiveEqnAllaire(
                                F_y_ptr,
                                vel_y_ptr,
                                Z_rho_B_ptr,
                                Z_rho_T_ptr,
                                m_B_ptr,
                                m_T_ptr,
                                &E[idx_cell_B],
                                &E[idx_cell_T],
                                Z_B_ptr,
                                Z_T_ptr,
                                Y_DIRECTION);
                            
                            // Mulitply fluxes by dt.
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                *F_y_ptr[ei] *= dt;
                            }
                        }
                    }
                    
                    // Compute the source.
                    // Get the grid spacing.
                    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                            patch.getPatchGeometry()));
                    
                    const double* dx = patch_geom->getDx();
                    
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        double *S = source->getPointer(d_num_eqn - (d_num_species - 1) + si);
                        
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0]; i++)
                            {
                                // Compute the indices of cell and faces. 
                                const int idx_cell_wghost = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_nghost = i + j*interior_dims[0];
                                
                                const int idx_face_x_L = i +
                                    j*(interior_dims[0] + 1);
                                
                                const int idx_face_x_R = (i + 1) +
                                    j*(interior_dims[0] + 1);
                                
                                const int idx_face_y_B = j +
                                    i*(interior_dims[1] + 1);
                                
                                const int idx_face_y_T = (j + 1) +
                                    i*(interior_dims[1] + 1);
                                
                                const int& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                                const int& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                                
                                const int& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                                const int& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                                
                                S[idx_cell_nghost] += dt*Z[si][idx_cell_wghost]*((u_R - u_L)/dx[0] + (v_T - v_B)/dx[1]);
                            }
                        }
                    }
                }  // if (d_dim == tbox::Dimension(2))
                else if (d_dim == tbox::Dimension(3))
                {
                    // Get the arrays of time-dependent variables.
                    std::vector<double*> Z_rho;
                    for (int si = 0; si < d_num_species; si++)
                    {
                        Z_rho.push_back(partial_density->getPointer(si));
                    }  
                    double* rho_u = momentum->getPointer(0);
                    double* rho_v = momentum->getPointer(1);
                    double* rho_w = momentum->getPointer(2);
                    double* E     = total_energy->getPointer(0);
                    std::vector<double*> Z;
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        Z.push_back(volume_fraction->getPointer(si));
                    }
                    
                    // Compute the fluxes in the x direction and velocity components at the face
                    // normal to the x direction.
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
                                
                                std::vector<double*> F_x_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_x_ptr.push_back(&(convective_flux->getPointer(0, ei)[idx_face_x]));
                                }
                                
                                std::vector<double*> vel_x_ptr;
                                for (int vi = 0; vi < 3; vi++)
                                {
                                    vel_x_ptr.push_back(&(velocity_intercell->getPointer(0, vi)[idx_face_x]));
                                }
                                
                                std::vector<const double*> Z_rho_L_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_L_ptr.push_back(
                                        &Z_rho[si][idx_cell_L]);
                                }
                                
                                std::vector<const double*> Z_rho_R_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_R_ptr.push_back(
                                        &Z_rho[si][idx_cell_R]);
                                }
                                
                                std::vector<const double*> m_L_ptr;
                                m_L_ptr.push_back(&rho_u[idx_cell_L]);
                                m_L_ptr.push_back(&rho_v[idx_cell_L]);
                                m_L_ptr.push_back(&rho_w[idx_cell_L]);
                                
                                std::vector<const double*> m_R_ptr;
                                m_R_ptr.push_back(&rho_u[idx_cell_R]);
                                m_R_ptr.push_back(&rho_v[idx_cell_R]);
                                m_R_ptr.push_back(&rho_w[idx_cell_R]);
                                
                                std::vector<const double*> Z_L_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_L_ptr.push_back(
                                        &Z[si][idx_cell_L]);
                                }
                                
                                std::vector<const double*> Z_R_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_R_ptr.push_back(
                                        &Z[si][idx_cell_R]);
                                }
                                
                                d_Riemann_solver.computeIntercellFluxAndVelocityForFiveEqnAllaire(
                                    F_x_ptr,
                                    vel_x_ptr,
                                    Z_rho_L_ptr,
                                    Z_rho_R_ptr,
                                    m_L_ptr,
                                    m_R_ptr,
                                    &E[idx_cell_L],
                                    &E[idx_cell_R],
                                    Z_L_ptr,
                                    Z_R_ptr,
                                    X_DIRECTION);
                                
                                // Mulitply fluxes by dt.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    *F_x_ptr[ei] *= dt;
                                }
                            }
                        }
                    }
                    
                    // Compute the fluxes in the y direction and velocity components at the face
                    // normal to the y direction.
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
                                
                                std::vector<double*> F_y_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_y_ptr.push_back(&(convective_flux->getPointer(1, ei)[idx_face_y]));
                                }
                                
                                std::vector<double*> vel_y_ptr;
                                for (int vi = 0; vi < 3; vi++)
                                {
                                    vel_y_ptr.push_back(&(velocity_intercell->getPointer(1, vi)[idx_face_y]));
                                }
                                
                                std::vector<const double*> Z_rho_B_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_B_ptr.push_back(
                                        &Z_rho[si][idx_cell_B]);
                                }
                                
                                std::vector<const double*> Z_rho_T_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_T_ptr.push_back(
                                        &Z_rho[si][idx_cell_T]);
                                }
                                
                                std::vector<const double*> m_B_ptr;
                                m_B_ptr.push_back(&rho_u[idx_cell_B]);
                                m_B_ptr.push_back(&rho_v[idx_cell_B]);
                                m_B_ptr.push_back(&rho_w[idx_cell_B]);
                                
                                std::vector<const double*> m_T_ptr;
                                m_T_ptr.push_back(&rho_u[idx_cell_T]);
                                m_T_ptr.push_back(&rho_v[idx_cell_T]);
                                m_T_ptr.push_back(&rho_w[idx_cell_T]);
                                
                                std::vector<const double*> Z_B_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_B_ptr.push_back(
                                        &Z[si][idx_cell_B]);
                                }
                                
                                std::vector<const double*> Z_T_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_T_ptr.push_back(
                                        &Z[si][idx_cell_T]);
                                }
                                
                                d_Riemann_solver.computeIntercellFluxAndVelocityForFiveEqnAllaire(
                                    F_y_ptr,
                                    vel_y_ptr,
                                    Z_rho_B_ptr,
                                    Z_rho_T_ptr,
                                    m_B_ptr,
                                    m_T_ptr,
                                    &E[idx_cell_B],
                                    &E[idx_cell_T],
                                    Z_B_ptr,
                                    Z_T_ptr,
                                    Y_DIRECTION);
                                
                                // Mulitply fluxes by dt.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    *F_y_ptr[ei] *= dt;
                                }
                            }
                        }
                    }
                    
                    // Compute the fluxes in the z direction and velocity components at the face
                    // normal to the z direction.
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
                                
                                std::vector<double*> F_z_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_z_ptr.push_back(&(convective_flux->getPointer(2, ei)[idx_face_z]));
                                }
                                
                                std::vector<double*> vel_z_ptr;
                                for (int vi = 0; vi < 3; vi++)
                                {
                                    vel_z_ptr.push_back(&(velocity_intercell->getPointer(2, vi)[idx_face_z]));
                                }
                                
                                std::vector<const double*> Z_rho_B_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_B_ptr.push_back(
                                        &Z_rho[si][idx_cell_B]);
                                }
                                
                                std::vector<const double*> Z_rho_F_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_F_ptr.push_back(
                                        &Z_rho[si][idx_cell_F]);
                                }
                                
                                std::vector<const double*> m_B_ptr;
                                m_B_ptr.push_back(&rho_u[idx_cell_B]);
                                m_B_ptr.push_back(&rho_v[idx_cell_B]);
                                m_B_ptr.push_back(&rho_w[idx_cell_B]);
                                
                                std::vector<const double*> m_F_ptr;
                                m_F_ptr.push_back(&rho_u[idx_cell_F]);
                                m_F_ptr.push_back(&rho_v[idx_cell_F]);
                                m_F_ptr.push_back(&rho_w[idx_cell_F]);
                                
                                std::vector<const double*> Z_B_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_B_ptr.push_back(
                                        &Z[si][idx_cell_B]);
                                }
                                
                                std::vector<const double*> Z_F_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_F_ptr.push_back(
                                        &Z[si][idx_cell_F]);
                                }
                                
                                d_Riemann_solver.computeIntercellFluxAndVelocityForFiveEqnAllaire(
                                    F_z_ptr,
                                    vel_z_ptr,
                                    Z_rho_B_ptr,
                                    Z_rho_F_ptr,
                                    m_B_ptr,
                                    m_F_ptr,
                                    &E[idx_cell_B],
                                    &E[idx_cell_F],
                                    Z_B_ptr,
                                    Z_F_ptr,
                                    Z_DIRECTION);
                                
                                // Mulitply fluxes by dt.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    *F_z_ptr[ei] *= dt;
                                }
                            }
                        }
                    }
                    
                    // Compute the source.
                    // Get the grid spacing.
                    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                            patch.getPatchGeometry()));
                    
                    const double* dx = patch_geom->getDx();
                    
                    for (int si = 0; si < d_num_species - 1; si++)
                    {
                        double *S = source->getPointer(d_num_eqn - (d_num_species - 1) + si);
                        
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = 0; j < interior_dims[1]; j++)
                            {
                                for (int i = 0; i < interior_dims[0]; i++)
                                {
                                    // Compute the indices of cell and faces. 
                                    const int idx_cell_wghost = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_nghost = i +
                                        j*interior_dims[0] +
                                        k*interior_dims[0]*interior_dims[1];
                                    
                                    const int idx_face_x_L = i +
                                        j*(interior_dims[0] + 1) +
                                        k*(interior_dims[0] + 1)*interior_dims[1];
                                    
                                    const int idx_face_x_R = (i + 1) +
                                        j*(interior_dims[0] + 1) +
                                        k*(interior_dims[0] + 1)*interior_dims[1];
                                    
                                    const int idx_face_y_B = j +
                                        k*(interior_dims[1] + 1) +
                                        i*(interior_dims[1] + 1)*interior_dims[2];
                                    
                                    const int idx_face_y_T = (j + 1) +
                                        k*(interior_dims[1] + 1) +
                                        i*(interior_dims[1] + 1)*interior_dims[2];
                                    
                                    const int idx_face_z_B = k +
                                        i*(interior_dims[2] + 1) +
                                        j*(interior_dims[2] + 1)*interior_dims[0];
                                    
                                    const int idx_face_z_F = (k + 1) +
                                        i*(interior_dims[2] + 1) +
                                        j*(interior_dims[2] + 1)*interior_dims[0];
                                    
                                    const int& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                                    const int& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                                    
                                    const int& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                                    const int& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                                    
                                    const int& w_B = velocity_intercell->getPointer(2, 2)[idx_face_z_B];
                                    const int& w_F = velocity_intercell->getPointer(2, 2)[idx_face_z_F];
                                    
                                    S[idx_cell_nghost] += dt*Z[si][idx_cell_wghost]*((u_R - u_L)/dx[0] +
                                        (v_T - v_B)/dx[1] + (w_F - w_B)/dx[2]);
                                }
                            }
                        }
                    }
                }  // if (d_dim == tbox::Dimension(3))
                
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
