#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorLLF.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

ConvectiveFluxReconstructorLLF::ConvectiveFluxReconstructorLLF(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geom,
    const FLOW_MODEL& flow_model,
    const int& num_eqn,
    const int& num_species,
    const boost::shared_ptr<EquationOfState>& equation_of_state,
    const boost::shared_ptr<tbox::Database>& shock_capturing_scheme_db):
        ConvectiveFluxReconstructor(object_name,
            dim,
            grid_geom,
            flow_model,
            num_eqn,
            num_species,
            equation_of_state,
            shock_capturing_scheme_db)
{
    d_num_conv_ghosts = hier::IntVector::getOne(d_dim);
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorLLF::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorLLF object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorLLF: this = "
       << (ConvectiveFluxReconstructorLLF *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_variables_set = "
       << d_variables_set
       << std::endl;
    os << "d_num_ghosts_set = "
       << d_num_ghosts_set
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorLLF::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putString("d_shock_capturing_scheme", "LLF");
}


/*
 * Compute the convective fluxes and sources due to hyperbolization
 * of the equations.
 */
void
ConvectiveFluxReconstructorLLF::computeConvectiveFluxesAndSources(
    hier::Patch& patch,
    const double time,
    const double dt,
    const int RK_step_number,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    NULL_USE(RK_step_number);
    
    if (d_variables_set == true)
    {
        if (d_num_ghosts_set == true)
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
                    
                    // Allocate temporary patch data.
                    boost::shared_ptr<pdat::CellData<double> > velocity(
                        new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > pressure(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > sound_speed(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > spectral_radius(
                        new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
                    
                    if (d_dim == tbox::Dimension(1))
                    {
                        // Get the arrays of time-dependent variables.
                        double* rho   = density->getPointer(0);
                        double* rho_u = momentum->getPointer(0);
                        double* E     = total_energy->getPointer(0);
                        
                        // Get the arrays of useful variables.
                        double* u     = velocity->getPointer(0);
                        double* p     = pressure->getPointer(0);
                        double* c     = sound_speed->getPointer(0);
                        double* sp_x  = spectral_radius->getPointer(0);
                        
                        // Create a vector of arrays to time-dependent variables.
                        std::vector<double*> Q;
                        Q.push_back(density->getPointer(0));
                        Q.push_back(momentum->getPointer(0));
                        Q.push_back(total_energy->getPointer(0));
                        
                        // Compute the field of pressure, sound speed and spectral radius.
                        for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                        {
                            // Compute index into linear data array.
                            const int idx = i + d_num_ghosts[0];
                            
                            u[idx] = rho_u[idx]/rho[idx];
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&rho_u[idx]);
                            
                            p[idx] = d_equation_of_state->getPressure(
                                &rho[idx],
                                m_ptr,
                                &E[idx]);
                            
                            c[idx] = d_equation_of_state->getSoundSpeedWithPressure(
                                &rho[idx],
                                &p[idx]);
                            
                            sp_x[idx] = fabs(u[idx]) + c[idx];
                        }
                        
                        // Compute the fluxes in the x direction.
                        for (int i = 0; i < interior_dims[0] + 1; i++)
                        {
                            // Compute the indices of left cell, right cell and face.
                            const int idx_cell_L = i - 1 + d_num_ghosts[0];
                            
                            const int idx_cell_R = i + d_num_ghosts[0];
                            
                            const int idx_face_x = i;
                            
                            const double alpha_x = fmax(sp_x[idx_cell_L], sp_x[idx_cell_R]);
                            
                            std::vector<double> F_x_L;
                            std::vector<double> F_x_R;
                            
                            F_x_L.push_back(rho_u[idx_cell_L]);
                            F_x_L.push_back(rho_u[idx_cell_L]*u[idx_cell_L] + p[idx_cell_L]);
                            F_x_L.push_back(u[idx_cell_L]*(E[idx_cell_L] + p[idx_cell_L]));
                            
                            F_x_R.push_back(rho_u[idx_cell_R]);
                            F_x_R.push_back(rho_u[idx_cell_R]*u[idx_cell_R] + p[idx_cell_R]);
                            F_x_R.push_back(u[idx_cell_R]*(E[idx_cell_R] + p[idx_cell_R]));
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                double* F_x = convective_flux->getPointer(0, ei);
                                F_x[idx_face_x] = 0.5*dt*(F_x_L[ei] + F_x_R[ei] - alpha_x*(Q[ei][idx_cell_R] - Q[ei][idx_cell_L]));
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
                        
                        // Get the arrays of useful variables.
                        double* u     = velocity->getPointer(0);
                        double* v     = velocity->getPointer(1);
                        double* p     = pressure->getPointer(0);
                        double* c     = sound_speed->getPointer(0);
                        double* sp_x  = spectral_radius->getPointer(0);
                        double* sp_y  = spectral_radius->getPointer(1);
                        
                        // Create a vector of arrays to time-dependent variables.
                        std::vector<double*> Q;
                        Q.push_back(density->getPointer(0));
                        Q.push_back(momentum->getPointer(0));
                        Q.push_back(momentum->getPointer(1));
                        Q.push_back(total_energy->getPointer(0));
                        
                        // Compute the field of pressure, sound speed and spectral radius.
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
                                
                                p[idx] = d_equation_of_state->getPressure(
                                    &rho[idx],
                                    m_ptr,
                                    &E[idx]);
                                
                                c[idx] = d_equation_of_state->getSoundSpeedWithPressure(
                                    &rho[idx],
                                    &p[idx]);
                                
                                sp_x[idx] = fabs(u[idx]) + c[idx];
                                sp_y[idx] = fabs(v[idx]) + c[idx];
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
                                
                                const double alpha_x = fmax(sp_x[idx_cell_L], sp_x[idx_cell_R]);
                                
                                std::vector<double> F_x_L;
                                std::vector<double> F_x_R;
                                
                                F_x_L.push_back(rho_u[idx_cell_L]);
                                F_x_L.push_back(rho_u[idx_cell_L]*u[idx_cell_L] + p[idx_cell_L]);
                                F_x_L.push_back(rho_u[idx_cell_L]*v[idx_cell_L]);
                                F_x_L.push_back(u[idx_cell_L]*(E[idx_cell_L] + p[idx_cell_L]));
                                
                                F_x_R.push_back(rho_u[idx_cell_R]);
                                F_x_R.push_back(rho_u[idx_cell_R]*u[idx_cell_R] + p[idx_cell_R]);
                                F_x_R.push_back(rho_u[idx_cell_R]*v[idx_cell_R]);
                                F_x_R.push_back(u[idx_cell_R]*(E[idx_cell_R] + p[idx_cell_R]));
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double* F_x = convective_flux->getPointer(0, ei);
                                    F_x[idx_face_x] = 0.5*dt*(F_x_L[ei] + F_x_R[ei] - alpha_x*(Q[ei][idx_cell_R] - Q[ei][idx_cell_L]));
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
                                
                                const double alpha_y = fmax(sp_y[idx_cell_B], sp_y[idx_cell_T]);
                                
                                std::vector<double> F_y_B;
                                std::vector<double> F_y_T;
                                
                                F_y_B.push_back(rho_v[idx_cell_B]);
                                F_y_B.push_back(rho_v[idx_cell_B]*u[idx_cell_B]);
                                F_y_B.push_back(rho_v[idx_cell_B]*v[idx_cell_B] + p[idx_cell_B]);
                                F_y_B.push_back(v[idx_cell_B]*(E[idx_cell_B] + p[idx_cell_B]));
                                
                                F_y_T.push_back(rho_v[idx_cell_T]);
                                F_y_T.push_back(rho_v[idx_cell_T]*u[idx_cell_T]);
                                F_y_T.push_back(rho_v[idx_cell_T]*v[idx_cell_T] + p[idx_cell_T]);
                                F_y_T.push_back(v[idx_cell_T]*(E[idx_cell_T] + p[idx_cell_T]));
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double* F_y = convective_flux->getPointer(1, ei);
                                    F_y[idx_face_y] = 0.5*dt*(F_y_B[ei] + F_y_T[ei] - alpha_y*(Q[ei][idx_cell_T] - Q[ei][idx_cell_B]));
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
                        
                        // Get the arrays of useful variables.
                        double* u     = velocity->getPointer(0);
                        double* v     = velocity->getPointer(1);
                        double* w     = velocity->getPointer(2);
                        double* p     = pressure->getPointer(0);
                        double* c     = sound_speed->getPointer(0);
                        double* sp_x  = spectral_radius->getPointer(0);
                        double* sp_y  = spectral_radius->getPointer(1);
                        double* sp_z  = spectral_radius->getPointer(2);
                        
                        // Create a vector of arrays to time-dependent variables.
                        std::vector<double*> Q;
                        Q.push_back(density->getPointer(0));
                        Q.push_back(momentum->getPointer(0));
                        Q.push_back(momentum->getPointer(1));
                        Q.push_back(momentum->getPointer(2));
                        Q.push_back(total_energy->getPointer(0));
                        
                        // Compute the field of pressure, sound speed and spectral radius.
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
                                    
                                    p[idx] = d_equation_of_state->getPressure(
                                        &rho[idx],
                                        m_ptr,
                                        &E[idx]);
                                    
                                    c[idx] = d_equation_of_state->getSoundSpeedWithPressure(
                                        &rho[idx],
                                        &p[idx]);
                                    
                                    sp_x[idx] = fabs(u[idx]) + c[idx];
                                    sp_y[idx] = fabs(v[idx]) + c[idx];
                                    sp_z[idx] = fabs(w[idx]) + c[idx];
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
                                    
                                    const double alpha_x = fmax(sp_x[idx_cell_L], sp_x[idx_cell_R]);
                                    
                                    std::vector<double> F_x_L;
                                    std::vector<double> F_x_R;
                                    
                                    F_x_L.push_back(rho_u[idx_cell_L]);
                                    F_x_L.push_back(rho_u[idx_cell_L]*u[idx_cell_L] + p[idx_cell_L]);
                                    F_x_L.push_back(rho_u[idx_cell_L]*v[idx_cell_L]);
                                    F_x_L.push_back(rho_u[idx_cell_L]*w[idx_cell_L]);
                                    F_x_L.push_back(u[idx_cell_L]*(E[idx_cell_L] + p[idx_cell_L]));
                                    
                                    F_x_R.push_back(rho_u[idx_cell_R]);
                                    F_x_R.push_back(rho_u[idx_cell_R]*u[idx_cell_R] + p[idx_cell_R]);
                                    F_x_R.push_back(rho_u[idx_cell_R]*v[idx_cell_R]);
                                    F_x_R.push_back(rho_u[idx_cell_R]*w[idx_cell_R]);
                                    F_x_R.push_back(u[idx_cell_R]*(E[idx_cell_R] + p[idx_cell_R]));
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_x = convective_flux->getPointer(0, ei);
                                        F_x[idx_face_x] = 0.5*dt*(F_x_L[ei] + F_x_R[ei] - alpha_x*(Q[ei][idx_cell_R] - Q[ei][idx_cell_L]));
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
                                    
                                    const double alpha_y = fmax(sp_y[idx_cell_B], sp_y[idx_cell_T]);
                                    
                                    std::vector<double> F_y_B;
                                    std::vector<double> F_y_T;
                                    
                                    F_y_B.push_back(rho_v[idx_cell_B]);
                                    F_y_B.push_back(rho_v[idx_cell_B]*u[idx_cell_B]);
                                    F_y_B.push_back(rho_v[idx_cell_B]*v[idx_cell_B] + p[idx_cell_B]);
                                    F_y_B.push_back(rho_v[idx_cell_B]*w[idx_cell_B]);
                                    F_y_B.push_back(v[idx_cell_B]*(E[idx_cell_B] + p[idx_cell_B]));
                                    
                                    F_y_T.push_back(rho_v[idx_cell_T]);
                                    F_y_T.push_back(rho_v[idx_cell_T]*u[idx_cell_T]);
                                    F_y_T.push_back(rho_v[idx_cell_T]*v[idx_cell_T] + p[idx_cell_T]);
                                    F_y_T.push_back(rho_v[idx_cell_T]*w[idx_cell_T]);
                                    F_y_T.push_back(v[idx_cell_T]*(E[idx_cell_T] + p[idx_cell_T]));
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_y = convective_flux->getPointer(1, ei);
                                        F_y[idx_face_y] = 0.5*dt*(F_y_B[ei] + F_y_T[ei] - alpha_y*(Q[ei][idx_cell_T] - Q[ei][idx_cell_B]));                    
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
                                    
                                    const double alpha_z = fmax(sp_z[idx_cell_B], sp_z[idx_cell_F]);
                                    
                                    std::vector<double> F_z_B;
                                    std::vector<double> F_z_F;
                                    
                                    F_z_B.push_back(rho_w[idx_cell_B]);
                                    F_z_B.push_back(rho_w[idx_cell_B]*u[idx_cell_B]);
                                    F_z_B.push_back(rho_w[idx_cell_B]*v[idx_cell_B]);
                                    F_z_B.push_back(rho_w[idx_cell_B]*w[idx_cell_B] + p[idx_cell_B]);
                                    F_z_B.push_back(w[idx_cell_B]*(E[idx_cell_B] + p[idx_cell_B]));
                                    
                                    F_z_F.push_back(rho_w[idx_cell_F]);
                                    F_z_F.push_back(rho_w[idx_cell_F]*u[idx_cell_F]);
                                    F_z_F.push_back(rho_w[idx_cell_F]*v[idx_cell_F]);
                                    F_z_F.push_back(rho_w[idx_cell_F]*w[idx_cell_F] + p[idx_cell_F]);
                                    F_z_F.push_back(w[idx_cell_F]*(E[idx_cell_F] + p[idx_cell_F]));
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_z = convective_flux->getPointer(2, ei);
                                        F_z[idx_face_z] = 0.5*dt*(F_z_B[ei] + F_z_F[ei] - alpha_z*(Q[ei][idx_cell_F] - Q[ei][idx_cell_B]));                    
                                    }
                                }
                            }
                        }
                    }   // if (d_dim == tbox::Dimension(3))
                    
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
                    
                    boost::shared_ptr<pdat::CellData<double> > velocity(
                        new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > pressure(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > sound_speed(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > spectral_radius(
                        new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
                    
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
                        
                        // Get the arrays of useful variables.
                        double* u     = velocity->getPointer(0);
                        double* p     = pressure->getPointer(0);
                        double* c     = sound_speed->getPointer(0);
                        double* sp_x  = spectral_radius->getPointer(0);
                        
                        // Create a vector of arrays to time-dependent variables.
                        std::vector<double*> Q;
                        Q.push_back(density->getPointer(0));
                        Q.push_back(momentum->getPointer(0));
                        Q.push_back(total_energy->getPointer(0));
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Q.push_back(mass_fraction->getPointer(si));
                        }
                        
                        // Compute the field of total density, pressure, sound speed and spectral radius.
                        for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                        {
                            // Compute index into linear data array.
                            const int idx = i + d_num_ghosts[0];
                            
                            std::vector<const double*> Y_idx;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Y_idx.push_back(&Y[si][idx]);
                            }
                            
                            u[idx] = rho_u[idx]/rho[idx];
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&rho_u[idx]);
                            
                            p[idx] = d_equation_of_state->
                                getPressureWithMassFraction(
                                    &rho[idx],
                                    m_ptr,
                                    &E[idx],
                                    Y_idx);
                            
                            c[idx] = d_equation_of_state->
                                getSoundSpeedWithMassFractionAndPressure(
                                    &rho[idx],
                                    Y_idx,
                                    &p[idx]);
                            
                            sp_x[idx] = fabs(u[idx]) + c[idx];
                        }
                        
                        // Compute the fluxes in the x direction and velocity components at the face
                        // normal to the x direction.
                        double* u_intercell = velocity_intercell->getPointer(0, 0);
                        for (int i = 0; i < interior_dims[0] + 1; i++)
                        {
                            // Compute the indices of left cell, right cell and face.
                            const int idx_cell_L = i - 1 + d_num_ghosts[0];
                            
                            const int idx_cell_R = i + d_num_ghosts[0];
                            
                            const int idx_face_x = i;
                            
                            const double alpha_x = fmax(sp_x[idx_cell_L], sp_x[idx_cell_R]);
                            
                            std::vector<double> F_x_L;
                            std::vector<double> F_x_R;
                            
                            F_x_L.push_back(rho_u[idx_cell_L]);
                            F_x_L.push_back(rho_u[idx_cell_L]*u[idx_cell_L] + p[idx_cell_L]);
                            F_x_L.push_back(u[idx_cell_L]*(E[idx_cell_L] + p[idx_cell_L]));
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_x_L.push_back(u[idx_cell_L]*Y[si][idx_cell_L]);
                            }
                            
                            F_x_R.push_back(rho_u[idx_cell_R]);
                            F_x_R.push_back(rho_u[idx_cell_R]*u[idx_cell_R] + p[idx_cell_R]);
                            F_x_R.push_back(u[idx_cell_R]*(E[idx_cell_R] + p[idx_cell_R]));
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_x_R.push_back(u[idx_cell_R]*Y[si][idx_cell_R]);
                            }
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                double *F_x = convective_flux->getPointer(0, ei);
                                F_x[idx_face_x] = 0.5*dt*(F_x_L[ei] + F_x_R[ei] -
                                    alpha_x*(Q[ei][idx_cell_R] - Q[ei][idx_cell_L]));
                            }
                            
                            u_intercell[idx_face_x] = (p[idx_cell_R] - p[idx_cell_L] +
                                rho_u[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                rho_u[idx_cell_R]*(alpha_x - u[idx_cell_R]))/
                                (rho[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                rho[idx_cell_R]*(alpha_x - u[idx_cell_R]));
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
                        
                        // Get the arrays of useful variables.
                        double* u     = velocity->getPointer(0);
                        double* v     = velocity->getPointer(1);
                        double* p     = pressure->getPointer(0);
                        double* c     = sound_speed->getPointer(0);
                        double* sp_x  = spectral_radius->getPointer(0);
                        double* sp_y  = spectral_radius->getPointer(1);
                        
                        // Create a vector of arrays to time-dependent variables.
                        std::vector<double*> Q;
                        Q.push_back(density->getPointer(0));
                        Q.push_back(momentum->getPointer(0));
                        Q.push_back(momentum->getPointer(1));
                        Q.push_back(total_energy->getPointer(0));
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Q.push_back(mass_fraction->getPointer(si));
                        }
                        
                        // Compute the field of total density, pressure, sound speed and spectral radius.
                        for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                        {
                            for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                            {
                                // Compute index into linear data array.
                                const int idx = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                std::vector<const double*> Y_idx;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_idx.push_back(&Y[si][idx]);
                                }
                                
                                u[idx] = rho_u[idx]/rho[idx];
                                v[idx] = rho_v[idx]/rho[idx];
                                
                                std::vector<const double*> m_ptr;
                                m_ptr.push_back(&rho_u[idx]);
                                m_ptr.push_back(&rho_v[idx]);
                                
                                p[idx] = d_equation_of_state->
                                    getPressureWithMassFraction(
                                        &rho[idx],
                                        m_ptr,
                                        &E[idx],
                                        Y_idx);
                                
                                c[idx] = d_equation_of_state->
                                    getSoundSpeedWithMassFractionAndPressure(
                                        &rho[idx],
                                        Y_idx,
                                        &p[idx]);
                                
                                sp_x[idx] = fabs(u[idx]) + c[idx];
                                sp_y[idx] = fabs(v[idx]) + c[idx];
                            }
                        }
                        
                        // Compute the fluxes in the x direction and velocity components at the face
                        // normal to the x direction.
                        double* u_intercell = velocity_intercell->getPointer(0, 0);
                        double* v_intercell = velocity_intercell->getPointer(0, 1);
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
                                
                                const double alpha_x = fmax(sp_x[idx_cell_L], sp_x[idx_cell_R]);
                                
                                std::vector<double> F_x_L;
                                std::vector<double> F_x_R;
                                
                                F_x_L.push_back(rho_u[idx_cell_L]);
                                F_x_L.push_back(rho_u[idx_cell_L]*u[idx_cell_L] + p[idx_cell_L]);
                                F_x_L.push_back(rho_u[idx_cell_L]*v[idx_cell_L]);
                                F_x_L.push_back(u[idx_cell_L]*(E[idx_cell_L] + p[idx_cell_L]));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_x_L.push_back(u[idx_cell_L]*Y[si][idx_cell_L]);
                                }
                                
                                F_x_R.push_back(rho_u[idx_cell_R]);
                                F_x_R.push_back(rho_u[idx_cell_R]*u[idx_cell_R] + p[idx_cell_R]);
                                F_x_R.push_back(rho_u[idx_cell_R]*v[idx_cell_R]);
                                F_x_R.push_back(u[idx_cell_R]*(E[idx_cell_R] + p[idx_cell_R]));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_x_R.push_back(u[idx_cell_R]*Y[si][idx_cell_R]);
                                }
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double *F_x = convective_flux->getPointer(0, ei);
                                    F_x[idx_face_x] = 0.5*dt*(F_x_L[ei] + F_x_R[ei] -
                                        alpha_x*(Q[ei][idx_cell_R] - Q[ei][idx_cell_L]));
                                }
                                
                                u_intercell[idx_face_x] = (p[idx_cell_R] - p[idx_cell_L] +
                                    rho_u[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                    rho_u[idx_cell_R]*(alpha_x - u[idx_cell_R]))/
                                    (rho[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                    rho[idx_cell_R]*(alpha_x - u[idx_cell_R]));
                                
                                v_intercell[idx_face_x] = (rho_v[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                    rho_v[idx_cell_R]*(alpha_x - u[idx_cell_R]))/
                                    (rho[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                    rho[idx_cell_R]*(alpha_x - u[idx_cell_R]));
                            }
                        }
                        
                        // Compute the fluxes in the y direction and velocity components at the face
                        // normal to the y direction.
                        u_intercell = velocity_intercell->getPointer(1, 0);
                        v_intercell = velocity_intercell->getPointer(1, 1);
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
                                
                                const double alpha_y = fmax(sp_y[idx_cell_B], sp_y[idx_cell_T]);
                                
                                std::vector<double> F_y_B;
                                std::vector<double> F_y_T;
                                
                                F_y_B.push_back(rho_v[idx_cell_B]);
                                F_y_B.push_back(rho_v[idx_cell_B]*u[idx_cell_B]);
                                F_y_B.push_back(rho_v[idx_cell_B]*v[idx_cell_B] + p[idx_cell_B]);
                                F_y_B.push_back(v[idx_cell_B]*(E[idx_cell_B] + p[idx_cell_B]));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_y_B.push_back(v[idx_cell_B]*Y[si][idx_cell_B]);
                                }
                                
                                F_y_T.push_back(rho_v[idx_cell_T]);
                                F_y_T.push_back(rho_v[idx_cell_T]*u[idx_cell_T]);
                                F_y_T.push_back(rho_v[idx_cell_T]*v[idx_cell_T] + p[idx_cell_T]);
                                F_y_T.push_back(v[idx_cell_T]*(E[idx_cell_T] + p[idx_cell_T]));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_y_T.push_back(v[idx_cell_T]*Y[si][idx_cell_T]);
                                }
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double *F_y = convective_flux->getPointer(1, ei);
                                    F_y[idx_face_y] = 0.5*dt*(F_y_B[ei] + F_y_T[ei] -
                                        alpha_y*(Q[ei][idx_cell_T] - Q[ei][idx_cell_B]));
                                }
                                
                                u_intercell[idx_face_y] = (rho_u[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                    rho_u[idx_cell_T]*(alpha_y - v[idx_cell_T]))/
                                    (rho[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                    rho[idx_cell_T]*(alpha_y - v[idx_cell_T]));
                                
                                v_intercell[idx_face_y] = (p[idx_cell_T] - p[idx_cell_B] +
                                    rho_v[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                    rho_v[idx_cell_T]*(alpha_y - v[idx_cell_T]))/
                                    (rho[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                    rho[idx_cell_T]*(alpha_y - v[idx_cell_T]));
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
                        
                        // Get the arrays of useful variables.
                        double* u     = velocity->getPointer(0);
                        double* v     = velocity->getPointer(1);
                        double* w     = velocity->getPointer(2);
                        double* p     = pressure->getPointer(0);
                        double* c     = sound_speed->getPointer(0);
                        double* sp_x  = spectral_radius->getPointer(0);
                        double* sp_y  = spectral_radius->getPointer(1);
                        double* sp_z  = spectral_radius->getPointer(2);
                        
                        // Create a vector of arrays to time-dependent variables.
                        std::vector<double*> Q;
                        Q.push_back(density->getPointer(0));
                        Q.push_back(momentum->getPointer(0));
                        Q.push_back(momentum->getPointer(1));
                        Q.push_back(momentum->getPointer(2));
                        Q.push_back(total_energy->getPointer(0));
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Q.push_back(mass_fraction->getPointer(si));
                        }
                        
                        // Compute the field of total density, pressure, sound speed and spectral radius.
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
                                    
                                    std::vector<const double*> Y_idx;
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        Y_idx.push_back(&Y[si][idx]);
                                    }
                                    
                                    u[idx] = rho_u[idx]/rho[idx];
                                    v[idx] = rho_v[idx]/rho[idx];
                                    w[idx] = rho_w[idx]/rho[idx];
                                    
                                    std::vector<const double*> m_ptr;
                                    m_ptr.push_back(&rho_u[idx]);
                                    m_ptr.push_back(&rho_v[idx]);
                                    m_ptr.push_back(&rho_w[idx]);
                                    
                                    p[idx] = d_equation_of_state->
                                        getPressureWithMassFraction(
                                            &rho[idx],
                                            m_ptr,
                                            &E[idx],
                                            Y_idx);
                                    
                                    c[idx] = d_equation_of_state->
                                        getSoundSpeedWithMassFractionAndPressure(
                                            &rho[idx],
                                            Y_idx,
                                            &p[idx]);
                                    
                                    sp_x[idx] = fabs(u[idx]) + c[idx];
                                    sp_y[idx] = fabs(v[idx]) + c[idx];
                                    sp_z[idx] = fabs(w[idx]) + c[idx];
                                }
                            }
                        }
                        
                        // Compute the fluxes in the x direction and velocity components at the face
                        // normal to the x direction.
                        double* u_intercell = velocity_intercell->getPointer(0, 0);
                        double* v_intercell = velocity_intercell->getPointer(0, 1);
                        double* w_intercell = velocity_intercell->getPointer(0, 2);
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
                                    
                                    const double alpha_x = fmax(sp_x[idx_cell_L], sp_x[idx_cell_R]);
                                    
                                    std::vector<double> F_x_L;
                                    std::vector<double> F_x_R;
                                    
                                    F_x_L.push_back(rho_u[idx_cell_L]);
                                    F_x_L.push_back(rho_u[idx_cell_L]*u[idx_cell_L] + p[idx_cell_L]);
                                    F_x_L.push_back(rho_u[idx_cell_L]*v[idx_cell_L]);
                                    F_x_L.push_back(rho_u[idx_cell_L]*w[idx_cell_L]);
                                    F_x_L.push_back(u[idx_cell_L]*(E[idx_cell_L] + p[idx_cell_L]));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_x_L.push_back(u[idx_cell_L]*Y[si][idx_cell_L]);
                                    }
                                    
                                    F_x_R.push_back(rho_u[idx_cell_R]);
                                    F_x_R.push_back(rho_u[idx_cell_R]*u[idx_cell_R] + p[idx_cell_R]);
                                    F_x_R.push_back(rho_u[idx_cell_R]*v[idx_cell_R]);
                                    F_x_R.push_back(rho_u[idx_cell_R]*w[idx_cell_R]);
                                    F_x_R.push_back(u[idx_cell_R]*(E[idx_cell_R] + p[idx_cell_R]));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_x_R.push_back(u[idx_cell_R]*Y[si][idx_cell_R]);
                                    }
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double *F_x = convective_flux->getPointer(0, ei);
                                        F_x[idx_face_x] = 0.5*dt*(F_x_L[ei] + F_x_R[ei] -
                                            alpha_x*(Q[ei][idx_cell_R] - Q[ei][idx_cell_L]));
                                    }
                                    
                                    u_intercell[idx_face_x] = (p[idx_cell_R] - p[idx_cell_L] +
                                        rho_u[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                        rho_u[idx_cell_R]*(alpha_x - u[idx_cell_R]))/
                                        (rho[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                        rho[idx_cell_R]*(alpha_x - u[idx_cell_R]));
                                    
                                    v_intercell[idx_face_x] = (rho_v[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                        rho_v[idx_cell_R]*(alpha_x - u[idx_cell_R]))/
                                        (rho[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                        rho[idx_cell_R]*(alpha_x - u[idx_cell_R]));
                                    
                                    w_intercell[idx_face_x] = (rho_w[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                        rho_w[idx_cell_R]*(alpha_x - u[idx_cell_R]))/
                                        (rho[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                        rho[idx_cell_R]*(alpha_x - u[idx_cell_R]));
                                }
                            }
                        }
                        
                        // Compute the fluxes in the y direction and velocity components at the face
                        // normal to the y direction.
                        u_intercell = velocity_intercell->getPointer(1, 0);
                        v_intercell = velocity_intercell->getPointer(1, 1);
                        w_intercell = velocity_intercell->getPointer(1, 2);
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
                                    
                                    const double alpha_y = fmax(sp_y[idx_cell_B], sp_y[idx_cell_T]);
                                    
                                    std::vector<double> F_y_B;
                                    std::vector<double> F_y_T;
                                    
                                    F_y_B.push_back(rho_v[idx_cell_B]);
                                    F_y_B.push_back(rho_v[idx_cell_B]*u[idx_cell_B]);
                                    F_y_B.push_back(rho_v[idx_cell_B]*v[idx_cell_B] + p[idx_cell_B]);
                                    F_y_B.push_back(rho_v[idx_cell_B]*w[idx_cell_B]);
                                    F_y_B.push_back(v[idx_cell_B]*(E[idx_cell_B] + p[idx_cell_B]));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_y_B.push_back(v[idx_cell_B]*Y[si][idx_cell_B]);
                                    }
                                    
                                    F_y_T.push_back(rho_v[idx_cell_T]);
                                    F_y_T.push_back(rho_v[idx_cell_T]*u[idx_cell_T]);
                                    F_y_T.push_back(rho_v[idx_cell_T]*v[idx_cell_T] + p[idx_cell_T]);
                                    F_y_T.push_back(rho_v[idx_cell_T]*w[idx_cell_T]);
                                    F_y_T.push_back(v[idx_cell_T]*(E[idx_cell_T] + p[idx_cell_T]));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_y_T.push_back(v[idx_cell_T]*Y[si][idx_cell_T]);
                                    }
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double *F_y = convective_flux->getPointer(1, ei);
                                        F_y[idx_face_y] = 0.5*dt*(F_y_B[ei] + F_y_T[ei] -
                                            alpha_y*(Q[ei][idx_cell_T] - Q[ei][idx_cell_B]));
                                    }
                                    
                                    u_intercell[idx_face_y] = (rho_u[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                        rho_u[idx_cell_T]*(alpha_y - v[idx_cell_T]))/
                                        (rho[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                        rho[idx_cell_T]*(alpha_y - v[idx_cell_T]));
                                    
                                    v_intercell[idx_face_y] = (p[idx_cell_T] - p[idx_cell_B] +
                                        rho_v[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                        rho_v[idx_cell_T]*(alpha_y - v[idx_cell_T]))/
                                        (rho[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                        rho[idx_cell_T]*(alpha_y - v[idx_cell_T]));
                                    
                                    w_intercell[idx_face_y] = (rho_w[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                        rho_w[idx_cell_T]*(alpha_y - v[idx_cell_T]))/
                                        (rho[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                        rho[idx_cell_T]*(alpha_y - v[idx_cell_T]));
                                }
                            }
                        }
                        
                        // Compute the fluxes in the z direction and velocity components at the face
                        // normal to the z direction.
                        u_intercell = velocity_intercell->getPointer(2, 0);
                        v_intercell = velocity_intercell->getPointer(2, 1);
                        w_intercell = velocity_intercell->getPointer(2, 2);
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
                                    
                                    const double alpha_z = fmax(sp_z[idx_cell_B], sp_z[idx_cell_F]);
                                    
                                    std::vector<double> F_z_B;
                                    std::vector<double> F_z_F;
                                    
                                    F_z_B.push_back(rho_w[idx_cell_B]);
                                    F_z_B.push_back(rho_w[idx_cell_B]*u[idx_cell_B]);
                                    F_z_B.push_back(rho_w[idx_cell_B]*v[idx_cell_B]);
                                    F_z_B.push_back(rho_w[idx_cell_B]*w[idx_cell_B] + p[idx_cell_B]);
                                    F_z_B.push_back(w[idx_cell_B]*(E[idx_cell_B] + p[idx_cell_B]));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_z_B.push_back(w[idx_cell_B]*Y[si][idx_cell_B]);
                                    }
                                    
                                    F_z_F.push_back(rho_w[idx_cell_F]);
                                    F_z_F.push_back(rho_w[idx_cell_F]*u[idx_cell_F]);
                                    F_z_F.push_back(rho_w[idx_cell_F]*v[idx_cell_F]);
                                    F_z_F.push_back(rho_w[idx_cell_F]*w[idx_cell_F] + p[idx_cell_F]);
                                    F_z_F.push_back(w[idx_cell_F]*(E[idx_cell_F] + p[idx_cell_F]));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_z_F.push_back(w[idx_cell_F]*Y[si][idx_cell_F]);
                                    }
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_z = convective_flux->getPointer(2, ei);
                                        F_z[idx_face_z] = 0.5*dt*(F_z_B[ei] + F_z_F[ei] -
                                            alpha_z*(Q[ei][idx_cell_F] - Q[ei][idx_cell_B]));                    
                                    }
                                    
                                    u_intercell[idx_face_z] = (rho_u[idx_cell_B]*(-alpha_z - w[idx_cell_B]) -
                                        rho_u[idx_cell_F]*(alpha_z - w[idx_cell_F]))/
                                        (rho[idx_cell_B]*(-alpha_z - w[idx_cell_B]) -
                                        rho[idx_cell_F]*(alpha_z - w[idx_cell_F]));
                                    
                                    v_intercell[idx_face_z] = (rho_v[idx_cell_B]*(-alpha_z - w[idx_cell_B]) -
                                        rho_v[idx_cell_F]*(alpha_z - w[idx_cell_F]))/
                                        (rho[idx_cell_B]*(-alpha_z - w[idx_cell_B]) -
                                        rho[idx_cell_F]*(alpha_z - w[idx_cell_F]));
                                    
                                    w_intercell[idx_face_z] = (p[idx_cell_F] - p[idx_cell_B] +
                                        rho_w[idx_cell_B]*(-alpha_z - w[idx_cell_B]) -
                                        rho_w[idx_cell_F]*(alpha_z - w[idx_cell_F]))/
                                        (rho[idx_cell_B]*(-alpha_z - w[idx_cell_B]) -
                                        rho[idx_cell_F]*(alpha_z - w[idx_cell_F]));
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
                    
                    boost::shared_ptr<pdat::CellData<double> > density(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > velocity(
                        new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > pressure(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > sound_speed(
                        new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                    
                    boost::shared_ptr<pdat::CellData<double> > spectral_radius(
                        new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
                    
                    
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
                        
                        // Get the arrays of useful variables.
                        double* rho   = density->getPointer(0);
                        double* u     = velocity->getPointer(0);
                        double* p     = pressure->getPointer(0);
                        double* c     = sound_speed->getPointer(0);
                        double* sp_x  = spectral_radius->getPointer(0);
                        
                        // Create a vector of arrays to time-dependent variables.
                        std::vector<double*> Q;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Q.push_back(partial_density->getPointer(si));
                        }
                        Q.push_back(momentum->getPointer(0));
                        Q.push_back(total_energy->getPointer(0));
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Q.push_back(volume_fraction->getPointer(si));
                        }
                        
                        // Compute the field of total density, pressure, sound speed and spectral radius.
                        for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                        {
                            // Compute index into linear data array.
                            const int idx = i + d_num_ghosts[0];
                            
                            std::vector<const double*> Z_rho_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_ptr.push_back(&Z_rho[si][idx]);
                            }
                            
                            std::vector<const double*> Z_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Z_ptr.push_back(&Z[si][idx]);
                            }
                            
                            rho[idx] = d_equation_of_state->
                                getTotalDensity(
                                    Z_rho_ptr);
                            
                            u[idx] = rho_u[idx]/rho[idx];
                            
                            std::vector<const double*> m_ptr;
                            m_ptr.push_back(&rho_u[idx]);
                            
                            p[idx] = d_equation_of_state->
                                getPressureWithVolumeFraction(
                                    &rho[idx],
                                    m_ptr,
                                    &E[idx],
                                    Z_ptr);
                            
                            c[idx] = d_equation_of_state->
                                getSoundSpeedWithVolumeFractionAndPressure(
                                    &rho[idx],
                                    Z_ptr,
                                    &p[idx]);
                            
                            sp_x[idx] = fabs(u[idx]) + c[idx];
                        }
                        
                        // Compute the fluxes in the x direction and velocity components at the face
                        // normal to the x direction.
                        double* u_intercell = velocity_intercell->getPointer(0, 0);
                        for (int i = 0; i < interior_dims[0] + 1; i++)
                        {
                            // Compute the indices of left cell, right cell and face.
                            const int idx_cell_L = i - 1 + d_num_ghosts[0];
                            
                            const int idx_cell_R = i + d_num_ghosts[0];
                            
                            const int idx_face_x = i;
                            
                            const double alpha_x = fmax(sp_x[idx_cell_L], sp_x[idx_cell_R]);
                            
                            std::vector<double> F_x_L;
                            std::vector<double> F_x_R;
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                F_x_L.push_back(u[idx_cell_L]*Z_rho[si][idx_cell_L]);
                            }
                            F_x_L.push_back(rho_u[idx_cell_L]*u[idx_cell_L] + p[idx_cell_L]);
                            F_x_L.push_back(u[idx_cell_L]*(E[idx_cell_L] + p[idx_cell_L]));
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_x_L.push_back(u[idx_cell_L]*Z[si][idx_cell_L]);
                            }
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                F_x_R.push_back(u[idx_cell_R]*Z_rho[si][idx_cell_R]);
                            }
                            F_x_R.push_back(rho_u[idx_cell_R]*u[idx_cell_R] + p[idx_cell_R]);
                            F_x_R.push_back(u[idx_cell_R]*(E[idx_cell_R] + p[idx_cell_R]));
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_x_R.push_back(u[idx_cell_R]*Z[si][idx_cell_R]);
                            }
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                double *F_x = convective_flux->getPointer(0, ei);
                                F_x[idx_face_x] = 0.5*dt*(F_x_L[ei] + F_x_R[ei] -
                                    alpha_x*(Q[ei][idx_cell_R] - Q[ei][idx_cell_L]));
                            }
                            
                            u_intercell[idx_face_x] = (p[idx_cell_R] - p[idx_cell_L] +
                                rho_u[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                rho_u[idx_cell_R]*(alpha_x - u[idx_cell_R]))/
                                (rho[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                rho[idx_cell_R]*(alpha_x - u[idx_cell_R]));
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
                    } // if (d_dim == tbox::Dimension(1))
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
                        
                        // Get the arrays of useful variables.
                        double* rho   = density->getPointer(0);
                        double* u     = velocity->getPointer(0);
                        double* v     = velocity->getPointer(1);
                        double* p     = pressure->getPointer(0);
                        double* c     = sound_speed->getPointer(0);
                        double* sp_x  = spectral_radius->getPointer(0);
                        double* sp_y  = spectral_radius->getPointer(1);
                        
                        // Create a vector of arrays to time-dependent variables.
                        std::vector<double*> Q;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Q.push_back(partial_density->getPointer(si));
                        }
                        Q.push_back(momentum->getPointer(0));
                        Q.push_back(momentum->getPointer(1));
                        Q.push_back(total_energy->getPointer(0));
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Q.push_back(volume_fraction->getPointer(si));
                        }
                        
                        // Compute the field of total density, pressure, sound speed and spectral radius.
                        for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                        {
                            for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                            {
                                // Compute index into linear data array.
                                const int idx = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                std::vector<const double*> Z_rho_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_ptr.push_back(&Z_rho[si][idx]);
                                }
                                
                                std::vector<const double*> Z_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_ptr.push_back(&Z[si][idx]);
                                }
                                
                                rho[idx] = d_equation_of_state->
                                    getTotalDensity(
                                        Z_rho_ptr);
                                
                                u[idx] = rho_u[idx]/rho[idx];
                                v[idx] = rho_v[idx]/rho[idx];
                                
                                std::vector<const double*> m_ptr;
                                m_ptr.push_back(&rho_u[idx]);
                                m_ptr.push_back(&rho_v[idx]);
                                
                                p[idx] = d_equation_of_state->
                                    getPressureWithVolumeFraction(
                                        &rho[idx],
                                        m_ptr,
                                        &E[idx],
                                        Z_ptr);
                                
                                c[idx] = d_equation_of_state->
                                    getSoundSpeedWithVolumeFractionAndPressure(
                                        &rho[idx],
                                        Z_ptr,
                                        &p[idx]);
                                
                                sp_x[idx] = fabs(u[idx]) + c[idx];
                                sp_y[idx] = fabs(v[idx]) + c[idx];
                            }
                        }
                        
                        // Compute the fluxes in the x direction and velocity components at the face
                        // normal to the x direction.
                        double* u_intercell = velocity_intercell->getPointer(0, 0);
                        double* v_intercell = velocity_intercell->getPointer(0, 1);
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
                                
                                const double alpha_x = fmax(sp_x[idx_cell_L], sp_x[idx_cell_R]);
                                
                                std::vector<double> F_x_L;
                                std::vector<double> F_x_R;
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_x_L.push_back(u[idx_cell_L]*Z_rho[si][idx_cell_L]);
                                }
                                F_x_L.push_back(rho_u[idx_cell_L]*u[idx_cell_L] + p[idx_cell_L]);
                                F_x_L.push_back(rho_u[idx_cell_L]*v[idx_cell_L]);
                                F_x_L.push_back(u[idx_cell_L]*(E[idx_cell_L] + p[idx_cell_L]));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_x_L.push_back(u[idx_cell_L]*Z[si][idx_cell_L]);
                                }
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_x_R.push_back(u[idx_cell_R]*Z_rho[si][idx_cell_R]);
                                }
                                F_x_R.push_back(rho_u[idx_cell_R]*u[idx_cell_R] + p[idx_cell_R]);
                                F_x_R.push_back(rho_u[idx_cell_R]*v[idx_cell_R]);
                                F_x_R.push_back(u[idx_cell_R]*(E[idx_cell_R] + p[idx_cell_R]));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_x_R.push_back(u[idx_cell_R]*Z[si][idx_cell_R]);
                                }
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double *F_x = convective_flux->getPointer(0, ei);
                                    F_x[idx_face_x] = 0.5*dt*(F_x_L[ei] + F_x_R[ei] -
                                        alpha_x*(Q[ei][idx_cell_R] - Q[ei][idx_cell_L]));
                                }
                                
                                u_intercell[idx_face_x] = (p[idx_cell_R] - p[idx_cell_L] +
                                    rho_u[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                    rho_u[idx_cell_R]*(alpha_x - u[idx_cell_R]))/
                                    (rho[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                    rho[idx_cell_R]*(alpha_x - u[idx_cell_R]));
                                
                                v_intercell[idx_face_x] = (rho_v[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                    rho_v[idx_cell_R]*(alpha_x - u[idx_cell_R]))/
                                    (rho[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                    rho[idx_cell_R]*(alpha_x - u[idx_cell_R]));
                            }
                        }
                        
                        // Compute the fluxes in the y direction and velocity components at the face
                        // normal to the y direction.
                        u_intercell = velocity_intercell->getPointer(1, 0);
                        v_intercell = velocity_intercell->getPointer(1, 1);
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int j = 0; j < interior_dims[1] + 1; j++)
                            {
                                // Compute the indices of bottom cell, top cell and fluxes.
                                const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_T = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_face_y = j +
                                    i*(interior_dims[1] + 1);
                                
                                const double alpha_y = fmax(sp_y[idx_cell_B], sp_y[idx_cell_T]);
                                
                                std::vector<double> F_y_B;
                                std::vector<double> F_y_T;
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_y_B.push_back(v[idx_cell_B]*Z_rho[si][idx_cell_B]);
                                }
                                F_y_B.push_back(rho_v[idx_cell_B]*u[idx_cell_B]);
                                F_y_B.push_back(rho_v[idx_cell_B]*v[idx_cell_B] + p[idx_cell_B]);
                                F_y_B.push_back(v[idx_cell_B]*(E[idx_cell_B] + p[idx_cell_B]));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_y_B.push_back(v[idx_cell_B]*Z[si][idx_cell_B]);
                                }
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_y_T.push_back(v[idx_cell_T]*Z_rho[si][idx_cell_T]);
                                }
                                F_y_T.push_back(rho_v[idx_cell_T]*u[idx_cell_T]);
                                F_y_T.push_back(rho_v[idx_cell_T]*v[idx_cell_T] + p[idx_cell_T]);
                                F_y_T.push_back(v[idx_cell_T]*(E[idx_cell_T] + p[idx_cell_T]));
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_y_T.push_back(v[idx_cell_T]*Z[si][idx_cell_T]);
                                }
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    double *F_y = convective_flux->getPointer(1, ei);
                                    F_y[idx_face_y] = 0.5*dt*(F_y_B[ei] + F_y_T[ei] -
                                        alpha_y*(Q[ei][idx_cell_T] - Q[ei][idx_cell_B]));
                                }
                                
                                u_intercell[idx_face_y] = (rho_u[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                    rho_u[idx_cell_T]*(alpha_y - v[idx_cell_T]))/
                                    (rho[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                    rho[idx_cell_T]*(alpha_y - v[idx_cell_T]));
                                
                                v_intercell[idx_face_y] = (p[idx_cell_T] - p[idx_cell_B] +
                                    rho_v[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                    rho_v[idx_cell_T]*(alpha_y - v[idx_cell_T]))/
                                    (rho[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                    rho[idx_cell_T]*(alpha_y - v[idx_cell_T]));
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
                                    
                                    S[idx_cell_nghost] += dt*Z[si][idx_cell_wghost]*((u_R - u_L)/dx[0] +
                                        (v_T - v_B)/dx[1]);
                                }
                            }
                        }
                    } // if (d_dim == tbox::Dimension(2))
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
                        
                        // Get the arrays of useful variables.
                        double* rho   = density->getPointer(0);
                        double* u     = velocity->getPointer(0);
                        double* v     = velocity->getPointer(1);
                        double* w     = velocity->getPointer(2);
                        double* p     = pressure->getPointer(0);
                        double* c     = sound_speed->getPointer(0);
                        double* sp_x  = spectral_radius->getPointer(0);
                        double* sp_y  = spectral_radius->getPointer(1);
                        double* sp_z  = spectral_radius->getPointer(2);
                        
                        // Create a vector of arrays to time-dependent variables.
                        std::vector<double*> Q;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Q.push_back(partial_density->getPointer(si));
                        }
                        Q.push_back(momentum->getPointer(0));
                        Q.push_back(momentum->getPointer(1));
                        Q.push_back(momentum->getPointer(2));
                        Q.push_back(total_energy->getPointer(0));
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Q.push_back(volume_fraction->getPointer(si));
                        }
                        
                        // Compute the field of total density, pressure, sound speed and spectral radius.
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
                                    
                                    std::vector<const double*> Z_rho_ptr;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho_ptr.push_back(&Z_rho[si][idx]);
                                    }
                                    
                                    std::vector<const double*> Z_ptr;
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        Z_ptr.push_back(&Z[si][idx]);
                                    }
                                    
                                    rho[idx] = d_equation_of_state->
                                        getTotalDensity(
                                            Z_rho_ptr);
                                    
                                    u[idx] = rho_u[idx]/rho[idx];
                                    v[idx] = rho_v[idx]/rho[idx];
                                    w[idx] = rho_w[idx]/rho[idx];
                                    
                                    std::vector<const double*> m_ptr;
                                    m_ptr.push_back(&rho_u[idx]);
                                    m_ptr.push_back(&rho_v[idx]);
                                    m_ptr.push_back(&rho_w[idx]);
                                    
                                    p[idx] = d_equation_of_state->
                                        getPressureWithVolumeFraction(
                                            &rho[idx],
                                            m_ptr,
                                            &E[idx],
                                            Z_ptr);
                                    
                                    c[idx] = d_equation_of_state->
                                        getSoundSpeedWithVolumeFractionAndPressure(
                                            &rho[idx],
                                            Z_ptr,
                                            &p[idx]);
                                    
                                    sp_x[idx] = fabs(u[idx]) + c[idx];
                                    sp_y[idx] = fabs(v[idx]) + c[idx];
                                    sp_z[idx] = fabs(w[idx]) + c[idx];
                                }
                            }
                        }
                        
                        // Compute the fluxes in the x direction and velocity components at the face
                        // normal to the x direction.
                        double* u_intercell = velocity_intercell->getPointer(0, 0);
                        double* v_intercell = velocity_intercell->getPointer(0, 1);
                        double* w_intercell = velocity_intercell->getPointer(0, 2);
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
                                    
                                    const double alpha_x = fmax(sp_x[idx_cell_L], sp_x[idx_cell_R]);
                                    
                                    std::vector<double> F_x_L;
                                    std::vector<double> F_x_R;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_x_L.push_back(u[idx_cell_L]*Z_rho[si][idx_cell_L]);
                                    }
                                    F_x_L.push_back(rho_u[idx_cell_L]*u[idx_cell_L] + p[idx_cell_L]);
                                    F_x_L.push_back(rho_u[idx_cell_L]*v[idx_cell_L]);
                                    F_x_L.push_back(rho_u[idx_cell_L]*w[idx_cell_L]);
                                    F_x_L.push_back(u[idx_cell_L]*(E[idx_cell_L] + p[idx_cell_L]));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_x_L.push_back(u[idx_cell_L]*Z[si][idx_cell_L]);
                                    }
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_x_R.push_back(u[idx_cell_R]*Z_rho[si][idx_cell_R]);
                                    }
                                    F_x_R.push_back(rho_u[idx_cell_R]*u[idx_cell_R] + p[idx_cell_R]);
                                    F_x_R.push_back(rho_u[idx_cell_R]*v[idx_cell_R]);
                                    F_x_R.push_back(rho_u[idx_cell_R]*w[idx_cell_R]);
                                    F_x_R.push_back(u[idx_cell_R]*(E[idx_cell_R] + p[idx_cell_R]));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_x_R.push_back(u[idx_cell_R]*Z[si][idx_cell_R]);
                                    }
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double *F_x = convective_flux->getPointer(0, ei);
                                        F_x[idx_face_x] = 0.5*dt*(F_x_L[ei] + F_x_R[ei] -
                                            alpha_x*(Q[ei][idx_cell_R] - Q[ei][idx_cell_L]));
                                    }
                                    
                                    u_intercell[idx_face_x] = (p[idx_cell_R] - p[idx_cell_L] +
                                        rho_u[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                        rho_u[idx_cell_R]*(alpha_x - u[idx_cell_R]))/
                                        (rho[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                        rho[idx_cell_R]*(alpha_x - u[idx_cell_R]));
                                    
                                    v_intercell[idx_face_x] = (rho_v[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                        rho_v[idx_cell_R]*(alpha_x - u[idx_cell_R]))/
                                        (rho[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                        rho[idx_cell_R]*(alpha_x - u[idx_cell_R]));
                                    
                                    w_intercell[idx_face_x] = (rho_w[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                        rho_w[idx_cell_R]*(alpha_x - u[idx_cell_R]))/
                                        (rho[idx_cell_L]*(-alpha_x - u[idx_cell_L]) -
                                        rho[idx_cell_R]*(alpha_x - u[idx_cell_R]));
                                }
                            }
                        }
                           
                        // Compute the fluxes in the y direction and velocity components at the face
                        // normal to the y direction.
                        u_intercell = velocity_intercell->getPointer(1, 0);
                        v_intercell = velocity_intercell->getPointer(1, 1);
                        w_intercell = velocity_intercell->getPointer(1, 2);
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
                                    
                                    const double alpha_y = fmax(sp_y[idx_cell_B], sp_y[idx_cell_T]);
                                    
                                    std::vector<double> F_y_B;
                                    std::vector<double> F_y_T;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_y_B.push_back(v[idx_cell_B]*Z_rho[si][idx_cell_B]);
                                    }
                                    F_y_B.push_back(rho_v[idx_cell_B]*u[idx_cell_B]);
                                    F_y_B.push_back(rho_v[idx_cell_B]*v[idx_cell_B] + p[idx_cell_B]);
                                    F_y_B.push_back(rho_v[idx_cell_B]*w[idx_cell_B]);
                                    F_y_B.push_back(v[idx_cell_B]*(E[idx_cell_B] + p[idx_cell_B]));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_y_B.push_back(v[idx_cell_B]*Z[si][idx_cell_B]);
                                    }
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_y_T.push_back(v[idx_cell_T]*Z_rho[si][idx_cell_T]);
                                    }
                                    F_y_T.push_back(rho_v[idx_cell_T]*u[idx_cell_T]);
                                    F_y_T.push_back(rho_v[idx_cell_T]*v[idx_cell_T] + p[idx_cell_T]);
                                    F_y_T.push_back(rho_v[idx_cell_T]*w[idx_cell_T]);
                                    F_y_T.push_back(v[idx_cell_T]*(E[idx_cell_T] + p[idx_cell_T]));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_y_T.push_back(v[idx_cell_T]*Z[si][idx_cell_T]);
                                    }
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double *F_y = convective_flux->getPointer(1, ei);
                                        F_y[idx_face_y] = 0.5*dt*(F_y_B[ei] + F_y_T[ei] -
                                            alpha_y*(Q[ei][idx_cell_T] - Q[ei][idx_cell_B]));
                                    }
                                    
                                    u_intercell[idx_face_y] = (rho_u[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                        rho_u[idx_cell_T]*(alpha_y - v[idx_cell_T]))/
                                        (rho[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                        rho[idx_cell_T]*(alpha_y - v[idx_cell_T]));
                                    
                                    v_intercell[idx_face_y] = (p[idx_cell_T] - p[idx_cell_B] +
                                        rho_v[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                        rho_v[idx_cell_T]*(alpha_y - v[idx_cell_T]))/
                                        (rho[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                        rho[idx_cell_T]*(alpha_y - v[idx_cell_T]));
                                    
                                    w_intercell[idx_face_y] = (rho_w[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                        rho_w[idx_cell_T]*(alpha_y - v[idx_cell_T]))/
                                        (rho[idx_cell_B]*(-alpha_y - v[idx_cell_B]) -
                                        rho[idx_cell_T]*(alpha_y - v[idx_cell_T]));
                                }
                            }
                        }
                        
                        // Compute the fluxes in the z direction and velocity components at the face
                        // normal to the z direction.
                        u_intercell = velocity_intercell->getPointer(2, 0);
                        v_intercell = velocity_intercell->getPointer(2, 1);
                        w_intercell = velocity_intercell->getPointer(2, 2);
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
                                    
                                    const double alpha_z = fmax(sp_z[idx_cell_B], sp_z[idx_cell_F]);
                                    
                                    std::vector<double> F_z_B;
                                    std::vector<double> F_z_F;
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_z_B.push_back(w[idx_cell_B]*Z_rho[si][idx_cell_B]);
                                    }
                                    F_z_B.push_back(rho_w[idx_cell_B]*u[idx_cell_B]);
                                    F_z_B.push_back(rho_w[idx_cell_B]*v[idx_cell_B]);
                                    F_z_B.push_back(rho_w[idx_cell_B]*w[idx_cell_B] + p[idx_cell_B]);
                                    F_z_B.push_back(w[idx_cell_B]*(E[idx_cell_B] + p[idx_cell_B]));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_z_B.push_back(w[idx_cell_B]*Z[si][idx_cell_B]);
                                    }
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        F_z_F.push_back(w[idx_cell_F]*Z_rho[si][idx_cell_F]);
                                    }
                                    F_z_F.push_back(rho_w[idx_cell_F]*u[idx_cell_F]);
                                    F_z_F.push_back(rho_w[idx_cell_F]*v[idx_cell_F]);
                                    F_z_F.push_back(rho_w[idx_cell_F]*w[idx_cell_F] + p[idx_cell_F]);
                                    F_z_F.push_back(w[idx_cell_F]*(E[idx_cell_F] + p[idx_cell_F]));
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        F_z_F.push_back(w[idx_cell_F]*Z[si][idx_cell_F]);
                                    }
                                    
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        double* F_z = convective_flux->getPointer(2, ei);
                                        F_z[idx_face_z] = 0.5*dt*(F_z_B[ei] + F_z_F[ei] -
                                            alpha_z*(Q[ei][idx_cell_F] - Q[ei][idx_cell_B]));                    
                                    }
                                    
                                    u_intercell[idx_face_z] = (rho_u[idx_cell_B]*(-alpha_z - w[idx_cell_B]) -
                                        rho_u[idx_cell_F]*(alpha_z - w[idx_cell_F]))/
                                        (rho[idx_cell_B]*(-alpha_z - w[idx_cell_B]) -
                                        rho[idx_cell_F]*(alpha_z - w[idx_cell_F]));
                                    
                                    v_intercell[idx_face_z] = (rho_v[idx_cell_B]*(-alpha_z - w[idx_cell_B]) -
                                        rho_v[idx_cell_F]*(alpha_z - w[idx_cell_F]))/
                                        (rho[idx_cell_B]*(-alpha_z - w[idx_cell_B]) -
                                        rho[idx_cell_F]*(alpha_z - w[idx_cell_F]));
                                    
                                    w_intercell[idx_face_z] = (p[idx_cell_F] - p[idx_cell_B] +
                                        rho_w[idx_cell_B]*(-alpha_z - w[idx_cell_B]) -
                                        rho_w[idx_cell_F]*(alpha_z - w[idx_cell_F]))/
                                        (rho[idx_cell_B]*(-alpha_z - w[idx_cell_B]) -
                                        rho[idx_cell_F]*(alpha_z - w[idx_cell_F]));
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
                    } // if (d_dim == tbox::Dimension(3))
                    
                    break;
                }
                default:
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "d_flow_model '"
                        << d_flow_model
                        << "' not yet implemented."
                        << std::endl);
                }
            }
        } // if (set_variables == true )
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Number of ghost cells is not set yet."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Variables are not set yet."
            << std::endl);
    }
}
