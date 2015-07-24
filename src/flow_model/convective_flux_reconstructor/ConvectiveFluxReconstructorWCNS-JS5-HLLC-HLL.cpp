#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWCNS-JS5-HLLC-HLL.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL::ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL(
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
            "HLLC-HLL Riemann solver",
            d_dim,
            d_num_eqn,
            d_num_species,
            d_equation_of_state)
{
    d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*4;
    
    d_weights_d.push_back(1.0/16);
    d_weights_d.push_back(10.0/16);
    d_weights_d.push_back(5.0/16);
    
    d_weights_c.resize(boost::extents[3][3]);
    d_weights_c[0][0] = 3.0/8;
    d_weights_c[0][1] = -5.0/4;
    d_weights_c[0][2] = 15.0/8;
    d_weights_c[1][0] = -1.0/8;
    d_weights_c[1][1] = 3.0/4;
    d_weights_c[1][2] = 3.0/8;
    d_weights_c[2][0] = 3.0/8;
    d_weights_c[2][1] = 3.0/4;
    d_weights_c[2][2] = -1.0/8;
}


/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL::printClassData(
    std::ostream& os) const
{
    os << "\nConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL::printClassData..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL: this = "
       << (ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
    os << "d_set_variables = "
       << d_set_variables
       << std::endl;
    
    os << std::endl;
    
    os << "End of ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL::printClassData"
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putString("d_shock_capturing_scheme", "WCNS_JS5_HLLC_HLL");
}


/*
 * Compute the convective fluxes and sources due to hyperbolization
 * of the equtions.
 */
void
ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL::computeConvectiveFluxAndSource(
    hier::Patch& patch,
    const double time,
    const double dt,
    const boost::shared_ptr<hier::VariableContext> data_context)
{
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
                
                // Allocate temporary patch data.
                boost::shared_ptr<pdat::CellData<double> > velocity(
                    new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
                
                boost::shared_ptr<pdat::CellData<double> > pressure(
                    new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                
                boost::shared_ptr<pdat::CellData<double> > sound_speed(
                    new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > convective_flux_node;
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    convective_flux_node.push_back(boost::make_shared<pdat::CellData<double> >(
                        interior_box, d_num_eqn, d_num_ghosts));
                }
                
                boost::shared_ptr<pdat::FaceData<double> > convective_flux_midpoint(
                    new pdat::FaceData<double>(interior_box, d_num_eqn, hier::IntVector::getOne(d_dim)));
                
                boost::shared_ptr<pdat::FaceData<double> > projection_matrix(
                    new pdat::FaceData<double>(
                        interior_box,
                        d_num_eqn*d_num_eqn,
                        hier::IntVector::getOne(d_dim)));
                
                boost::shared_ptr<pdat::FaceData<double> > projection_matrix_inv(
                    new pdat::FaceData<double>(
                        interior_box,
                        d_num_eqn*d_num_eqn,
                        hier::IntVector::getOne(d_dim)));
                
                
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
                    std::vector<double*> F_x_node;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
                    }
                    std::vector<double*> F_x_midpoint;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_midpoint.push_back(convective_flux_midpoint->getPointer(0, ei));
                    }
                    
                    // Compute the field of velocities and pressure.
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
                        
                        F_x_node[0][idx] = rho_u[idx];
                        F_x_node[1][idx] = rho_u[idx]*u[idx] + p[idx];
                        F_x_node[3][idx] = u[idx]*(E[idx] + p[idx]);
                    }
                    
                    /*
                     * Compute the projection matrix and its inverse at the face normal to the
                     * x direction.
                     */
                    for (int i = -1; i < interior_dims[0] + 2; i++)
                    {
                        // Compute the indices of left cell, right cell and face of
                        // projection matrix.
                        const int idx_cell_L = i - 1 + d_num_ghosts[0];
                        
                        const int idx_cell_R = i + d_num_ghosts[0];
                        
                        const int idx_face_x = i + 1;
                        
                        // Get left and right quantities.
                        const double& rho_L = rho[idx_cell_L];
                        const double& rho_R = rho[idx_cell_R];
                        
                        const double& c_L = c[idx_cell_L];
                        const double& c_R = c[idx_cell_R];
                        
                        // Compute simply-averaged quantities.
                        const double rho_average = 0.5*(rho_L + rho_R);
                        const double c_average = 0.5*(c_L + c_R);
                        
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
                        
                        *R_x_intercell[0][0] = 1.0/(c_average*c_average);
                        *R_x_intercell[0][1] = 1.0;
                        *R_x_intercell[0][2] = 1.0/(c_average*c_average);
                        *R_x_intercell[1][0] = -1.0/(rho_average*c_average);
                        *R_x_intercell[1][1] = 0.0;
                        *R_x_intercell[1][2] = 1.0/(rho_average*c_average);
                        *R_x_intercell[2][0] = 1.0;
                        *R_x_intercell[2][1] = 0.0;
                        *R_x_intercell[2][2] = 1.0;
                        
                        *R_x_inv_intercell[0][0] = 0.0;
                        *R_x_inv_intercell[0][1] = -0.5*rho_average*c_average;
                        *R_x_inv_intercell[0][2] = 0.5;
                        *R_x_inv_intercell[1][0] = 1.0;
                        *R_x_inv_intercell[1][1] = 0.0;
                        *R_x_inv_intercell[1][2] = -1.0/(c_average*c_average);
                        *R_x_inv_intercell[2][0] = 0.0;
                        *R_x_inv_intercell[2][1] = 0.5*rho_average*c_average;
                        *R_x_inv_intercell[2][2] = 0.5;
                    }
                    
                    // Compute the mid-point flux in the x-direction.
                    for (int i = -1; i < interior_dims[0] + 2; i++)
                    {
                        // Compute the index of face of mid-point flux and
                        // projection matrix.
                        const int idx_face_x = i + 1;
                        
                        boost::multi_array<double, 2> W_array(
                            boost::extents[6][d_num_eqn]);
                        
                        /*
                         * Project primitive variables onto characteristic fields.
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
                            
                            std::vector<const double*> V_ptr;
                            V_ptr.push_back(&rho[idx_cell]);
                            V_ptr.push_back(&u[idx_cell]);
                            V_ptr.push_back(&p[idx_cell]);
                            
                            std::vector<double*> W_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                W_ptr.push_back(&W_array[m][ei]);
                            }
                            projectPrimitiveVariablesToCharacteristicFields(
                                W_ptr,
                                V_ptr,
                                R_x_inv_intercell);
                        }
                        
                        /*
                         * Do WENO interplation on characteristic variables to get W_L
                         * and W_R
                         */
                        
                        std::vector<double> W_L;
                        std::vector<double> W_R;
                        
                        performWENOInterpolation(W_L, W_R, W_array);
                        
                        /*
                         * Project characteristic variables back to physical fields.
                         */
                        
                        double rho_L;
                        double rho_R;
                        
                        std::vector<double> vel_L;
                        std::vector<double> vel_R;
                        vel_L.resize(d_dim.getValue());
                        vel_R.resize(d_dim.getValue());
                        
                        double p_L;
                        double p_R;
                        
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
                        
                        std::vector<const double*> W_L_ptr;
                        std::vector<const double*> W_R_ptr;
                        
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            W_L_ptr.push_back(&W_L[ei]);
                            W_R_ptr.push_back(&W_R[ei]);
                        }
                        
                        std::vector<double*> V_L_ptr;
                        V_L_ptr.push_back(&rho_L);
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            V_L_ptr.push_back(&vel_L[di]);
                        }
                        V_L_ptr.push_back(&p_L);
                        
                        std::vector<double*> V_R_ptr;
                        V_R_ptr.push_back(&rho_R);
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            V_R_ptr.push_back(&vel_R[di]);
                        }
                        V_R_ptr.push_back(&p_R);
                        
                        projectCharacteristicVariablesToPhysicalFields(
                            V_L_ptr,
                            W_L_ptr,
                            R_x_intercell);
                        
                        projectCharacteristicVariablesToPhysicalFields(
                            V_R_ptr,
                            W_R_ptr,
                            R_x_intercell);
                        
                        /*
                         * Convert the primitive variables into conservative variables.
                         */
                        
                        std::vector<double> m_L;
                        std::vector<double> m_R;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            m_L.push_back(rho_L*vel_L[di]);
                            m_R.push_back(rho_R*vel_R[di]);
                        }
                        
                        std::vector<const double*> vel_L_ptr;
                        std::vector<const double*> vel_R_ptr;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            vel_L_ptr.push_back(&vel_L[di]);
                            vel_R_ptr.push_back(&vel_R[di]);
                        }
                        
                        double E_L = d_equation_of_state->
                            getTotalEnergy(
                                &rho_L,
                                vel_L_ptr,
                                &p_L);
                        
                        double E_R = d_equation_of_state->
                            getTotalEnergy(
                                &rho_R,
                                vel_R_ptr,
                                &p_R);
                        
                        /*
                         * If the WENO interpolated density or energy are negative,
                         * use constant interpolation.
                         */
                        if ((rho_L < 0) || (rho_R < 0) || (E_L < 0) || (E_R < 0))
                        {
                            // Compute the indices of left cell and right cell.
                            const int idx_cell_L = i - 1 + d_num_ghosts[0];
                            const int idx_cell_R = i + d_num_ghosts[0];
                            
                            rho_L = rho[idx_cell_L];
                            rho_R = rho[idx_cell_R];
                            
                            m_L[0] = rho_u[idx_cell_L];
                            m_R[0] = rho_u[idx_cell_R];
                            
                            E_L = E[idx_cell_L];
                            E_R = E[idx_cell_R];
                        }
                        
                        /*
                         * Apply the Riemann solver.
                         */
                        
                        std::vector<const double*> m_L_ptr;
                        std::vector<const double*> m_R_ptr;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            m_L_ptr.push_back(&m_L[di]);
                            m_R_ptr.push_back(&m_R[di]);
                        }
                        
                        std::vector<double*> F_x_midpoint_ptr;
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            F_x_midpoint_ptr.push_back(&F_x_midpoint[ei][idx_face_x]);
                        }
                        
                        d_Riemann_solver.computeIntercellFluxForSingleSpecies(
                            F_x_midpoint_ptr,
                            &rho_L,
                            &rho_R,
                            m_L_ptr,
                            m_R_ptr,
                            &E_L,
                            &E_R,
                            X_DIRECTION);
                    }
                    
                    // Compute the flux in the x direction.
                    for (int i = 0; i < interior_dims[0] + 1; i++)
                    {
                        // Compute the indices.
                        const int idx_face_x = i;
                        const int idx_midpoint_x = i + 1;
                        const int idx_node_L = i - 1 + d_num_ghosts[0];
                        const int idx_node_R = i + d_num_ghosts[0];
                        
                        // Compute the flux.
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            convective_flux->getPointer(0, ei)[idx_face_x] = dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x + 1] +
                                F_x_midpoint[ei][idx_midpoint_x - 1]) -
                                3.0/10*(F_x_node[ei][idx_node_R] +
                                F_x_node[ei][idx_node_L]) +
                                23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
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
                    
                    // Get the arrays of temporary patch data.
                    double* u     = velocity->getPointer(0);
                    double* v     = velocity->getPointer(1);
                    double* p     = pressure->getPointer(0);
                    double* c     = sound_speed->getPointer(0);
                    std::vector<double*> F_x_node;
                    std::vector<double*> F_y_node;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
                        F_y_node.push_back(convective_flux_node[1]->getPointer(ei));
                    }
                    std::vector<double*> F_x_midpoint;
                    std::vector<double*> F_y_midpoint;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_midpoint.push_back(convective_flux_midpoint->getPointer(0, ei));
                        F_y_midpoint.push_back(convective_flux_midpoint->getPointer(1, ei));
                    }
                    
                    // Compute the field of velocities and pressure.
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
                            
                            F_x_node[0][idx] = rho_u[idx];
                            F_x_node[1][idx] = rho_u[idx]*u[idx] + p[idx];
                            F_x_node[2][idx] = rho_u[idx]*v[idx];
                            F_x_node[3][idx] = u[idx]*(E[idx] + p[idx]);
                            
                            F_y_node[0][idx] = rho_v[idx];
                            F_y_node[1][idx] = rho_v[idx]*u[idx];
                            F_y_node[2][idx] = rho_v[idx]*v[idx] + p[idx];
                            F_y_node[3][idx] = v[idx]*(E[idx] + p[idx]);
                        }
                    }
                    
                    /*
                     * Compute the projection matrix and its inverse at the face normal to the
                     * x direction.
                     */
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = -1; i < interior_dims[0] + 2; i++)
                        {
                            // Compute the indices of left cell, right cell and face of
                            // projection matrix.
                            const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_cell_R = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_face_x = (i + 1) +
                                (j + 1)*(interior_dims[0] + 3);
                            
                            // Get left and right quantities.
                            const double& rho_L = rho[idx_cell_L];
                            const double& rho_R = rho[idx_cell_R];
                            
                            const double& c_L = c[idx_cell_L];
                            const double& c_R = c[idx_cell_R];
                            
                            // Compute simply-averaged quantities.
                            const double rho_average = 0.5*(rho_L + rho_R);
                            const double c_average = 0.5*(c_L + c_R);
                            
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
                            
                            *R_x_intercell[0][0] = 1.0/(c_average*c_average);
                            *R_x_intercell[0][1] = 1.0;
                            *R_x_intercell[0][2] = 0.0;
                            *R_x_intercell[0][3] = 1.0/(c_average*c_average);
                            *R_x_intercell[1][0] = -1.0/(rho_average*c_average);
                            *R_x_intercell[1][1] = 0.0;
                            *R_x_intercell[1][2] = 0.0;
                            *R_x_intercell[1][3] = 1.0/(rho_average*c_average);
                            *R_x_intercell[2][0] = 0.0;
                            *R_x_intercell[2][1] = 0.0;
                            *R_x_intercell[2][2] = 1.0;
                            *R_x_intercell[2][3] = 0.0;
                            *R_x_intercell[3][0] = 1.0;
                            *R_x_intercell[3][1] = 0.0;
                            *R_x_intercell[3][2] = 0.0;
                            *R_x_intercell[3][3] = 1.0;
                            
                            *R_x_inv_intercell[0][0] = 0.0;
                            *R_x_inv_intercell[0][1] = -0.5*rho_average*c_average;
                            *R_x_inv_intercell[0][2] = 0.0;
                            *R_x_inv_intercell[0][3] = 0.5;
                            *R_x_inv_intercell[1][0] = 1.0;
                            *R_x_inv_intercell[1][1] = 0.0;
                            *R_x_inv_intercell[1][2] = 0.0;
                            *R_x_inv_intercell[1][3] = -1.0/(c_average*c_average);
                            *R_x_inv_intercell[2][0] = 0.0;
                            *R_x_inv_intercell[2][1] = 0.0;
                            *R_x_inv_intercell[2][2] = 1.0;
                            *R_x_inv_intercell[2][3] = 0.0;
                            *R_x_inv_intercell[3][0] = 0.0;
                            *R_x_inv_intercell[3][1] = 0.5*rho_average*c_average;
                            *R_x_inv_intercell[3][2] = 0.0;
                            *R_x_inv_intercell[3][3] = 0.5;
                        }
                    }
                    
                    /*
                     * Compute the projection matrix and its inverse at the face normal to the
                     * y direction.
                     */
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int j = -1; j < interior_dims[1] + 2; j++)
                        {
                            // Compute the indices of bottom cell, top cell and face of
                            // projection matrix.
                            const int idx_cell_B = (i + d_num_ghosts[0]) +
                                (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_cell_T = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_face_y = (j + 1) +
                                (i + 1)*(interior_dims[1] + 3);
                            
                            // Get bottom and top quantities.
                            
                            const double& rho_B = rho[idx_cell_B];
                            const double& rho_T = rho[idx_cell_T];
                            
                            const double& c_B = c[idx_cell_B];
                            const double& c_T = c[idx_cell_T];
                            
                            // Compute simply-averaged quantities.
                            const double rho_average = 0.5*(rho_B + rho_T);
                            const double c_average = 0.5*(c_B + c_T);
                            
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
                            
                            *R_y_intercell[0][0] = 1.0/(c_average*c_average);
                            *R_y_intercell[0][1] = 1.0;
                            *R_y_intercell[0][2] = 0.0;
                            *R_y_intercell[0][3] = 1.0/(c_average*c_average);
                            *R_y_intercell[1][0] = 0.0;
                            *R_y_intercell[1][1] = 0.0;
                            *R_y_intercell[1][2] = 1.0;
                            *R_y_intercell[1][3] = 0.0;
                            *R_y_intercell[2][0] = -1.0/(rho_average*c_average);
                            *R_y_intercell[2][1] = 0.0;
                            *R_y_intercell[2][2] = 0.0;
                            *R_y_intercell[2][3] = 1.0/(rho_average*c_average);
                            *R_y_intercell[3][0] = 1.0;
                            *R_y_intercell[3][1] = 0.0;
                            *R_y_intercell[3][2] = 0.0;
                            *R_y_intercell[3][3] = 1.0;
                            
                            *R_y_inv_intercell[0][0] = 0.0;
                            *R_y_inv_intercell[0][1] = 0.0;
                            *R_y_inv_intercell[0][2] = -0.5*rho_average*c_average;
                            *R_y_inv_intercell[0][3] = 0.5;
                            *R_y_inv_intercell[1][0] = 1.0;
                            *R_y_inv_intercell[1][1] = 0.0;
                            *R_y_inv_intercell[1][2] = 0.0;
                            *R_y_inv_intercell[1][3] = -1.0/(c_average*c_average);
                            *R_y_inv_intercell[2][0] = 0.0;
                            *R_y_inv_intercell[2][1] = 1.0;
                            *R_y_inv_intercell[2][2] = 0.0;
                            *R_y_inv_intercell[2][3] = 0.0;
                            *R_y_inv_intercell[3][0] = 0.0;
                            *R_y_inv_intercell[3][1] = 0.0;
                            *R_y_inv_intercell[3][2] = 0.5*rho_average*c_average;
                            *R_y_inv_intercell[3][3] = 0.5;
                        }
                    }
                    
                    // Compute the mid-point flux in the x direction.
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = -1; i < interior_dims[0] + 2; i++)
                        {
                            // Compute the index of face of mid-point flux and
                            // projection matrix.
                            const int idx_face_x = (i + 1) +
                                (j + 1)*(interior_dims[0] + 3);
                            
                            boost::multi_array<double, 2> W_array(
                                boost::extents[6][d_num_eqn]);
                            
                            /*
                             * Project primitive variables onto characteristic fields.
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
                                
                                std::vector<const double*> V_ptr;
                                V_ptr.push_back(&rho[idx_cell]);
                                V_ptr.push_back(&u[idx_cell]);
                                V_ptr.push_back(&v[idx_cell]);
                                V_ptr.push_back(&p[idx_cell]);
                                
                                std::vector<double*> W_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_ptr.push_back(&W_array[m][ei]);
                                }
                                projectPrimitiveVariablesToCharacteristicFields(
                                    W_ptr,
                                    V_ptr,
                                    R_x_inv_intercell);
                            }
                            
                            /*
                             * Do WENO interplation on characteristic variables to get W_L
                             * and W_R
                             */
                            
                            std::vector<double> W_L;
                            std::vector<double> W_R;
                            
                            performWENOInterpolation(W_L, W_R, W_array);
                            
                            /*
                             * Project characteristic variables back to physical fields.
                             */
                            
                            double rho_L;
                            double rho_R;
                            
                            std::vector<double> vel_L;
                            std::vector<double> vel_R;
                            vel_L.resize(d_dim.getValue());
                            vel_R.resize(d_dim.getValue());
                            
                            double p_L;
                            double p_R;
                            
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
                            
                            std::vector<const double*> W_L_ptr;
                            std::vector<const double*> W_R_ptr;
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                W_L_ptr.push_back(&W_L[ei]);
                                W_R_ptr.push_back(&W_R[ei]);
                            }
                            
                            std::vector<double*> V_L_ptr;
                            V_L_ptr.push_back(&rho_L);
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                V_L_ptr.push_back(&vel_L[di]);
                            }
                            V_L_ptr.push_back(&p_L);
                            
                            std::vector<double*> V_R_ptr;
                            V_R_ptr.push_back(&rho_R);
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                V_R_ptr.push_back(&vel_R[di]);
                            }
                            V_R_ptr.push_back(&p_R);
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                V_L_ptr,
                                W_L_ptr,
                                R_x_intercell);
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                V_R_ptr,
                                W_R_ptr,
                                R_x_intercell);
                            
                            /*
                             * Convert the primitive variables into conservative variables.
                             */
                            
                            std::vector<double> m_L;
                            std::vector<double> m_R;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_L.push_back(rho_L*vel_L[di]);
                                m_R.push_back(rho_R*vel_R[di]);
                            }
                            
                            std::vector<const double*> vel_L_ptr;
                            std::vector<const double*> vel_R_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                vel_L_ptr.push_back(&vel_L[di]);
                                vel_R_ptr.push_back(&vel_R[di]);
                            }
                            
                            double E_L = d_equation_of_state->
                                getTotalEnergy(
                                    &rho_L,
                                    vel_L_ptr,
                                    &p_L);
                            
                            double E_R = d_equation_of_state->
                                getTotalEnergy(
                                    &rho_R,
                                    vel_R_ptr,
                                    &p_R);
                            
                            /*
                             * If the WENO interpolated density or energy are negative,
                             * use constant interpolation.
                             */
                            if ((rho_L < 0) || (rho_R < 0) || (E_L < 0) || (E_R < 0))
                            {
                                // Compute the indices of left cell and right cell.
                                const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_R = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                rho_L = rho[idx_cell_L];
                                rho_R = rho[idx_cell_R];
                                
                                m_L[0] = rho_u[idx_cell_L];
                                m_L[1] = rho_v[idx_cell_L];
                                m_R[0] = rho_u[idx_cell_R];
                                m_R[1] = rho_v[idx_cell_R];
                                
                                E_L = E[idx_cell_L];
                                E_R = E[idx_cell_R];
                            }
                            
                            /*
                             * Apply the Riemann solver.
                             */
                            
                            std::vector<const double*> m_L_ptr;
                            std::vector<const double*> m_R_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_L_ptr.push_back(&m_L[di]);
                                m_R_ptr.push_back(&m_R[di]);
                            }
                            
                            std::vector<double*> F_x_midpoint_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_x_midpoint_ptr.push_back(&F_x_midpoint[ei][idx_face_x]);
                            }
                            
                            d_Riemann_solver.computeIntercellFluxForSingleSpecies(
                                F_x_midpoint_ptr,
                                &rho_L,
                                &rho_R,
                                m_L_ptr,
                                m_R_ptr,
                                &E_L,
                                &E_R,
                                X_DIRECTION);
                        }
                    }
                    
                    // Compute the mid-point flux in the y direction.
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int j = -1; j < interior_dims[1] + 2; j++)
                        {
                            // Compute the index of face of mid-point flux and
                            // projection matrix.
                            const int idx_face_y = (j + 1) +
                                (i + 1)*(interior_dims[1] + 3);
                            
                            boost::multi_array<double, 2> W_array(
                                boost::extents[6][d_num_eqn]);
                            
                            /*
                             * Project primitive variables onto characteristic fields.
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
                                
                                std::vector<const double*> V_ptr;
                                V_ptr.push_back(&rho[idx_cell]);
                                V_ptr.push_back(&u[idx_cell]);
                                V_ptr.push_back(&v[idx_cell]);
                                V_ptr.push_back(&p[idx_cell]);
                                
                                std::vector<double*> W_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_ptr.push_back(&W_array[m][ei]);
                                }
                                
                                projectPrimitiveVariablesToCharacteristicFields(
                                    W_ptr,
                                    V_ptr,
                                    R_y_inv_intercell);
                            }
                            
                            /*
                             * Do WENO interplation on characteristic variables to get W_B
                             * and W_T
                             */
                            
                            std::vector<double> W_B;
                            std::vector<double> W_T;
                            
                            performWENOInterpolation(W_B, W_T, W_array);
                            
                            /*
                             * Project characteristic variables back to physical fields.
                             */
                            
                            double rho_B;
                            double rho_T;
                            
                            std::vector<double> vel_B;
                            std::vector<double> vel_T;
                            vel_B.resize(d_dim.getValue());
                            vel_T.resize(d_dim.getValue());
                            
                            double p_B;
                            double p_T;
                            
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
                            
                            std::vector<const double*> W_B_ptr;
                            std::vector<const double*> W_T_ptr;
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                W_B_ptr.push_back(&W_B[ei]);
                                W_T_ptr.push_back(&W_T[ei]);
                            }
                            
                            std::vector<double*> V_B_ptr;
                            V_B_ptr.push_back(&rho_B);
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                V_B_ptr.push_back(&vel_B[di]);
                            }
                            V_B_ptr.push_back(&p_B);
                            
                            std::vector<double*> V_T_ptr;
                            V_T_ptr.push_back(&rho_T);
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                V_T_ptr.push_back(&vel_T[di]);
                            }
                            V_T_ptr.push_back(&p_T);
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                V_B_ptr,
                                W_B_ptr,
                                R_y_intercell);
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                V_T_ptr,
                                W_T_ptr,
                                R_y_intercell);
                            
                            /*
                             * Convert the primitive variables into conservative variables.
                             */
                            
                            std::vector<double> m_B;
                            std::vector<double> m_T;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_B.push_back(rho_B*vel_B[di]);
                                m_T.push_back(rho_T*vel_T[di]);
                            }
                            
                            std::vector<const double*> vel_B_ptr;
                            std::vector<const double*> vel_T_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                vel_B_ptr.push_back(&vel_B[di]);
                                vel_T_ptr.push_back(&vel_T[di]);
                            }
                            
                            double E_B = d_equation_of_state->
                                getTotalEnergy(
                                    &rho_B,
                                    vel_B_ptr,
                                    &p_B);
                            
                            double E_T = d_equation_of_state->
                                getTotalEnergy(
                                    &rho_T,
                                    vel_T_ptr,
                                    &p_T);
                            
                            /*
                             * If the WENO interpolated density or energy are negative,
                             * use constant interpolation.
                             */
                            if ((rho_B < 0) || (rho_T < 0) || (E_B < 0) || (E_T < 0))
                            {
                                // Compute the indices of bottom cell and top cell.
                                const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_T = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                rho_B = rho[idx_cell_B];
                                rho_T = rho[idx_cell_T];
                                
                                m_B[0] = rho_u[idx_cell_B];
                                m_B[1] = rho_v[idx_cell_B];
                                m_T[0] = rho_u[idx_cell_T];
                                m_T[1] = rho_v[idx_cell_T];
                                
                                E_B = E[idx_cell_B];
                                E_T = E[idx_cell_T];
                            }
                            
                            /*
                             * Apply the Riemann solver.
                             */
                            
                            std::vector<const double*> m_B_ptr;
                            std::vector<const double*> m_T_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_B_ptr.push_back(&m_B[di]);
                                m_T_ptr.push_back(&m_T[di]);
                            }
                            
                            std::vector<double*> F_y_midpoint_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_y_midpoint_ptr.push_back(&F_y_midpoint[ei][idx_face_y]);
                            }
                            
                            d_Riemann_solver.computeIntercellFluxForSingleSpecies(
                                F_y_midpoint_ptr,
                                &rho_B,
                                &rho_T,
                                m_B_ptr,
                                m_T_ptr,
                                &E_B,
                                &E_T,
                                Y_DIRECTION);
                        }
                    }
                    
                    // Compute the flux in the x direction.
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0] + 1; i++)
                        {
                            // Compute the indices.
                            const int idx_face_x = i +
                                j*(interior_dims[0] + 1);
                            
                            const int idx_midpoint_x = (i + 1) +
                                (j + 1)*(interior_dims[0] + 3);
                            
                            const int idx_node_L = (i - 1 + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_node_R = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            // Compute the flux.
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                convective_flux->getPointer(0, ei)[idx_face_x] = dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x + 1] +
                                    F_x_midpoint[ei][idx_midpoint_x - 1]) -
                                    3.0/10*(F_x_node[ei][idx_node_R] +
                                    F_x_node[ei][idx_node_L]) +
                                    23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
                            }
                        }
                    }
                    
                    // Compute the flux in the y direction.
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int j = 0; j < interior_dims[1] + 1; j++)
                        {
                            // Compute the indices.
                            const int idx_face_y = j +
                                i*(interior_dims[1] + 1);
                            
                            const int idx_midpoint_y = (j + 1) +
                                (i + 1)*(interior_dims[1] + 3);
                            
                            const int idx_node_B = (i + d_num_ghosts[0]) +
                                (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_node_T = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            // Compute the flux.
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                convective_flux->getPointer(1, ei)[idx_face_y] =
                                    dt*(1.0/30*(F_y_midpoint[ei][idx_midpoint_y + 1] +
                                    F_y_midpoint[ei][idx_midpoint_y - 1]) -
                                    3.0/10*(F_y_node[ei][idx_node_T] +
                                    F_y_node[ei][idx_node_B]) +
                                    23.0/15*F_y_midpoint[ei][idx_midpoint_y]);
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
                    
                    // Get the arrays of temporary patch data.
                    double* u     = velocity->getPointer(0);
                    double* v     = velocity->getPointer(1);
                    double* w     = velocity->getPointer(2);
                    double* p     = pressure->getPointer(0);
                    double* c     = sound_speed->getPointer(0);
                    std::vector<double*> F_x_node;
                    std::vector<double*> F_y_node;
                    std::vector<double*> F_z_node;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
                        F_y_node.push_back(convective_flux_node[1]->getPointer(ei));
                        F_z_node.push_back(convective_flux_node[2]->getPointer(ei));
                    }
                    std::vector<double*> F_x_midpoint;
                    std::vector<double*> F_y_midpoint;
                    std::vector<double*> F_z_midpoint;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_midpoint.push_back(convective_flux_midpoint->getPointer(0, ei));
                        F_y_midpoint.push_back(convective_flux_midpoint->getPointer(1, ei));
                        F_z_midpoint.push_back(convective_flux_midpoint->getPointer(2, ei));
                    }
                    
                    // Compute the field of velocities and pressure.
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
                            for (int i = -1; i < interior_dims[0] + 2; i++)
                            {
                                // Compute the indices of left cell, right cell and face of
                                // projection matrix.
                                const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_cell_R = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_face_x = (i + 1) +
                                    (j + 1)*(interior_dims[0] + 3) +
                                    (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                
                                // Get left and right quantities.
                                const double& rho_L = rho[idx_cell_L];
                                const double& rho_R = rho[idx_cell_R];
                                
                                const double& c_L = c[idx_cell_L];
                                const double& c_R = c[idx_cell_R];
                                
                                // Compute simply-averaged quantities.
                                const double rho_average = 0.5*(rho_L + rho_R);
                                const double c_average = 0.5*(c_L + c_R);
                                
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
                                
                                *R_x_intercell[0][0] = 1.0/(c_average*c_average);
                                *R_x_intercell[0][1] = 1.0;
                                *R_x_intercell[0][2] = 0.0;
                                *R_x_intercell[0][3] = 0.0;
                                *R_x_intercell[0][4] = 1.0/(c_average*c_average);
                                *R_x_intercell[1][0] = -1.0/(rho_average*c_average);
                                *R_x_intercell[1][1] = 0.0;
                                *R_x_intercell[1][2] = 0.0;
                                *R_x_intercell[1][3] = 0.0;
                                *R_x_intercell[1][4] = 1.0/(rho_average*c_average);
                                *R_x_intercell[2][0] = 0.0;
                                *R_x_intercell[2][1] = 0.0;
                                *R_x_intercell[2][2] = 1.0;
                                *R_x_intercell[2][3] = 0.0;
                                *R_x_intercell[2][4] = 0.0;
                                *R_x_intercell[3][0] = 0.0;
                                *R_x_intercell[3][1] = 0.0;
                                *R_x_intercell[3][2] = 0.0;
                                *R_x_intercell[3][3] = 1.0;
                                *R_x_intercell[3][4] = 0.0;
                                *R_x_intercell[4][0] = 1.0;
                                *R_x_intercell[4][1] = 0.0;
                                *R_x_intercell[4][2] = 0.0;
                                *R_x_intercell[4][3] = 0.0;
                                *R_x_intercell[4][4] = 1.0;
                                
                                *R_x_inv_intercell[0][0] = 0.0;
                                *R_x_inv_intercell[0][1] = -0.5*rho_average*c_average;
                                *R_x_inv_intercell[0][2] = 0.0;
                                *R_x_inv_intercell[0][3] = 0.0;
                                *R_x_inv_intercell[0][4] = 0.5;
                                *R_x_inv_intercell[1][0] = 1.0;
                                *R_x_inv_intercell[1][1] = 0.0;
                                *R_x_inv_intercell[1][2] = 0.0;
                                *R_x_inv_intercell[1][3] = 0.0;
                                *R_x_inv_intercell[1][4] = -1.0/(c_average*c_average);
                                *R_x_inv_intercell[2][0] = 0.0;
                                *R_x_inv_intercell[2][1] = 0.0;
                                *R_x_inv_intercell[2][2] = 1.0;
                                *R_x_inv_intercell[2][3] = 0.0;
                                *R_x_inv_intercell[2][4] = 0.0;
                                *R_x_inv_intercell[3][0] = 0.0;
                                *R_x_inv_intercell[3][1] = 0.0;
                                *R_x_inv_intercell[3][2] = 0.0;
                                *R_x_inv_intercell[3][3] = 1.0;
                                *R_x_inv_intercell[3][4] = 0.0;
                                *R_x_inv_intercell[4][0] = 0.0;
                                *R_x_inv_intercell[4][1] = 0.5*rho_average*c_average;
                                *R_x_inv_intercell[4][2] = 0.0;
                                *R_x_inv_intercell[4][3] = 0.0;
                                *R_x_inv_intercell[4][4] = 0.5;
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
                            for (int j = -1; j < interior_dims[1] + 2; j++)
                            {
                                // Compute the indices of bottom cell, top cell and face of
                                // projection matrix.
                                const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_cell_T = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_face_y = (j + 1) +
                                    (k + 1)*(interior_dims[1] + 3) +
                                    (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                
                                // Get bottom and top quantities.
                                
                                const double& rho_B = rho[idx_cell_B];
                                const double& rho_T = rho[idx_cell_T];
                                
                                const double& c_B = c[idx_cell_B];
                                const double& c_T = c[idx_cell_T];
                                
                                // Compute simply-averaged quantities.
                                const double rho_average = 0.5*(rho_B + rho_T);
                                const double c_average = 0.5*(c_B + c_T);
                                
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
                                
                                *R_y_intercell[0][0] = 1.0/(c_average*c_average);
                                *R_y_intercell[0][1] = 1.0;
                                *R_y_intercell[0][2] = 0.0;
                                *R_y_intercell[0][3] = 0.0;
                                *R_y_intercell[0][4] = 1.0/(c_average*c_average);
                                *R_y_intercell[1][0] = 0.0;
                                *R_y_intercell[1][1] = 0.0;
                                *R_y_intercell[1][2] = 1.0;
                                *R_y_intercell[1][3] = 0.0;
                                *R_y_intercell[1][4] = 0.0;
                                *R_y_intercell[2][0] = -1.0/(rho_average*c_average);
                                *R_y_intercell[2][1] = 0.0;
                                *R_y_intercell[2][2] = 0.0;
                                *R_y_intercell[2][3] = 0.0;
                                *R_y_intercell[2][4] = 1.0/(rho_average*c_average);
                                *R_y_intercell[3][0] = 0.0;
                                *R_y_intercell[3][1] = 0.0;
                                *R_y_intercell[3][2] = 0.0;
                                *R_y_intercell[3][3] = 1.0;
                                *R_y_intercell[3][4] = 0.0;
                                *R_y_intercell[4][0] = 1.0;
                                *R_y_intercell[4][1] = 0.0;
                                *R_y_intercell[4][2] = 0.0;
                                *R_y_intercell[4][3] = 0.0;
                                *R_y_intercell[4][4] = 1.0;
                                
                                *R_y_inv_intercell[0][0] = 0.0;
                                *R_y_inv_intercell[0][1] = 0.0;
                                *R_y_inv_intercell[0][2] = -0.5*rho_average*c_average;
                                *R_y_inv_intercell[0][3] = 0.0;
                                *R_y_inv_intercell[0][4] = 0.5;
                                *R_y_inv_intercell[1][0] = 1.0;
                                *R_y_inv_intercell[1][1] = 0.0;
                                *R_y_inv_intercell[1][2] = 0.0;
                                *R_y_inv_intercell[1][3] = 0.0;
                                *R_y_inv_intercell[1][4] = -1.0/(c_average*c_average);
                                *R_y_inv_intercell[2][0] = 0.0;
                                *R_y_inv_intercell[2][1] = 1.0;
                                *R_y_inv_intercell[2][2] = 0.0;
                                *R_y_inv_intercell[2][3] = 0.0;
                                *R_y_inv_intercell[2][4] = 0.0;
                                *R_y_inv_intercell[3][0] = 0.0;
                                *R_y_inv_intercell[3][1] = 0.0;
                                *R_y_inv_intercell[3][2] = 0.0;
                                *R_y_inv_intercell[3][3] = 1.0;
                                *R_y_inv_intercell[3][4] = 0.0;
                                *R_y_inv_intercell[4][0] = 0.0;
                                *R_y_inv_intercell[4][1] = 0.0;
                                *R_y_inv_intercell[4][2] = 0.5*rho_average*c_average;
                                *R_y_inv_intercell[4][3] = 0.0;
                                *R_y_inv_intercell[4][4] = 0.5;
                            }
                        }
                    }
                    
                    /*
                     * Compute the projection matrix and its inverse at the face normal to the
                     * z direction.
                     */
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int k = -1; k < interior_dims[2] + 2; k++)
                            {
                                // Compute the indices of back cell, front cell and face of
                                // projection matrix.
                                const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_cell_F = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_face_z = (k + 1) +
                                    (i + 1)*(interior_dims[2] + 3) +
                                    (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                
                                // Get back and front quantities.
                                
                                const double& rho_B = rho[idx_cell_B];
                                const double& rho_F = rho[idx_cell_F];
                                
                                const double& c_B = c[idx_cell_B];
                                const double& c_F = c[idx_cell_F];
                                
                                // Compute simply-averaged quantities.
                                const double rho_average = 0.5*(rho_B + rho_F);
                                const double c_average = 0.5*(c_B + c_F);
                                
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
                                
                                *R_z_intercell[0][0] = 1.0/(c_average*c_average);
                                *R_z_intercell[0][1] = 1.0;
                                *R_z_intercell[0][2] = 0.0;
                                *R_z_intercell[0][3] = 0.0;
                                *R_z_intercell[0][4] = 1.0/(c_average*c_average);
                                *R_z_intercell[1][0] = 0.0;
                                *R_z_intercell[1][1] = 0.0;
                                *R_z_intercell[1][2] = 1.0;
                                *R_z_intercell[1][3] = 0.0;
                                *R_z_intercell[1][4] = 0.0;
                                *R_z_intercell[2][0] = 0.0;
                                *R_z_intercell[2][1] = 0.0;
                                *R_z_intercell[2][2] = 0.0;
                                *R_z_intercell[2][3] = 1.0;
                                *R_z_intercell[2][4] = 0.0;
                                *R_z_intercell[3][0] = -1.0/(rho_average*c_average);
                                *R_z_intercell[3][1] = 0.0;
                                *R_z_intercell[3][2] = 0.0;
                                *R_z_intercell[3][3] = 0.0;
                                *R_z_intercell[3][4] = 1.0/(rho_average*c_average);
                                *R_z_intercell[4][0] = 1.0;
                                *R_z_intercell[4][1] = 0.0;
                                *R_z_intercell[4][2] = 0.0;
                                *R_z_intercell[4][3] = 0.0;
                                *R_z_intercell[4][4] = 1.0;
                                
                                *R_z_inv_intercell[0][0] = 0.0;
                                *R_z_inv_intercell[0][1] = 0.0;
                                *R_z_inv_intercell[0][2] = 0.0;
                                *R_z_inv_intercell[0][3] = -0.5*rho_average*c_average;
                                *R_z_inv_intercell[0][4] = 0.5;
                                *R_z_inv_intercell[1][0] = 1.0;
                                *R_z_inv_intercell[1][1] = 0.0;
                                *R_z_inv_intercell[1][2] = 0.0;
                                *R_z_inv_intercell[1][3] = 0.0;
                                *R_z_inv_intercell[1][4] = -1.0/(c_average*c_average);
                                *R_z_inv_intercell[2][0] = 0.0;
                                *R_z_inv_intercell[2][1] = 1.0;
                                *R_z_inv_intercell[2][2] = 0.0;
                                *R_z_inv_intercell[2][3] = 0.0;
                                *R_z_inv_intercell[2][4] = 0.0;
                                *R_z_inv_intercell[3][0] = 0.0;
                                *R_z_inv_intercell[3][1] = 0.0;
                                *R_z_inv_intercell[3][2] = 1.0;
                                *R_z_inv_intercell[3][3] = 0.0;
                                *R_z_inv_intercell[3][4] = 0.0;
                                *R_z_inv_intercell[4][0] = 0.0;
                                *R_z_inv_intercell[4][1] = 0.0;
                                *R_z_inv_intercell[4][2] = 0.0;
                                *R_z_inv_intercell[4][3] = 0.5*rho_average*c_average;
                                *R_z_inv_intercell[4][4] = 0.5;
                            }
                        }
                    }
                    
                    // Compute the mid-point flux in the x direction.
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = -1; i < interior_dims[0] + 2; i++)
                            {
                                // Compute the index of face of mid-point flux and
                                // projection matrix.
                                const int idx_face_x = (i + 1) +
                                    (j + 1)*(interior_dims[0] + 3) +
                                    (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                
                                boost::multi_array<double, 2> W_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                /*
                                 * Project primitive variables onto characteristic fields.
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
                                    
                                    std::vector<const double*> V_ptr;
                                    V_ptr.push_back(&rho[idx_cell]);
                                    V_ptr.push_back(&u[idx_cell]);
                                    V_ptr.push_back(&v[idx_cell]);
                                    V_ptr.push_back(&w[idx_cell]);
                                    V_ptr.push_back(&p[idx_cell]);
                                    
                                    std::vector<double*> W_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        W_ptr.push_back(&W_array[m][ei]);
                                    }
                                    projectPrimitiveVariablesToCharacteristicFields(
                                        W_ptr,
                                        V_ptr,
                                        R_x_inv_intercell);
                                }
                                
                                /*
                                 * Do WENO interplation on characteristic variables to get W_L
                                 * and W_R
                                 */
                                
                                std::vector<double> W_L;
                                std::vector<double> W_R;
                                
                                performWENOInterpolation(W_L, W_R, W_array);
                                
                                /*
                                 * Project characteristic variables back to physical fields.
                                 */
                                
                                double rho_L;
                                double rho_R;
                                
                                std::vector<double> vel_L;
                                std::vector<double> vel_R;
                                vel_L.resize(d_dim.getValue());
                                vel_R.resize(d_dim.getValue());
                                
                                double p_L;
                                double p_R;
                                
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
                                
                                std::vector<const double*> W_L_ptr;
                                std::vector<const double*> W_R_ptr;
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_L_ptr.push_back(&W_L[ei]);
                                    W_R_ptr.push_back(&W_R[ei]);
                                }
                                
                                std::vector<double*> V_L_ptr;
                                V_L_ptr.push_back(&rho_L);
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_L_ptr.push_back(&vel_L[di]);
                                }
                                V_L_ptr.push_back(&p_L);
                                
                                std::vector<double*> V_R_ptr;
                                V_R_ptr.push_back(&rho_R);
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_R_ptr.push_back(&vel_R[di]);
                                }
                                V_R_ptr.push_back(&p_R);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_L_ptr,
                                    W_L_ptr,
                                    R_x_intercell);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_R_ptr,
                                    W_R_ptr,
                                    R_x_intercell);
                                
                                /*
                                 * Convert the primitive variables into conservative variables.
                                 */
                                
                                std::vector<double> m_L;
                                std::vector<double> m_R;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_L.push_back(rho_L*vel_L[di]);
                                    m_R.push_back(rho_R*vel_R[di]);
                                }
                                
                                std::vector<const double*> vel_L_ptr;
                                std::vector<const double*> vel_R_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    vel_L_ptr.push_back(&vel_L[di]);
                                    vel_R_ptr.push_back(&vel_R[di]);
                                }
                                
                                double E_L = d_equation_of_state->
                                    getTotalEnergy(
                                        &rho_L,
                                        vel_L_ptr,
                                        &p_L);
                                
                                double E_R = d_equation_of_state->
                                    getTotalEnergy(
                                        &rho_R,
                                        vel_R_ptr,
                                        &p_R);
                                
                                /*
                                 * If the WENO interpolated density or energy are negative,
                                 * use constant interpolation.
                                 */
                                if ((rho_L < 0) || (rho_R < 0) || (E_L < 0) || (E_R < 0))
                                {
                                    // Compute the indices of left cell and right cell.
                                    const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_R = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    rho_L = rho[idx_cell_L];
                                    rho_R = rho[idx_cell_R];
                                    
                                    m_L[0] = rho_u[idx_cell_L];
                                    m_L[1] = rho_v[idx_cell_L];
                                    m_L[2] = rho_w[idx_cell_L];
                                    m_R[0] = rho_u[idx_cell_R];
                                    m_R[1] = rho_v[idx_cell_R];
                                    m_R[2] = rho_w[idx_cell_R];
                                    
                                    E_L = E[idx_cell_L];
                                    E_R = E[idx_cell_R];
                                }
                                
                                /*
                                 * Apply the Riemann solver.
                                 */
                                
                                std::vector<const double*> m_L_ptr;
                                std::vector<const double*> m_R_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_L_ptr.push_back(&m_L[di]);
                                    m_R_ptr.push_back(&m_R[di]);
                                }
                                
                                std::vector<double*> F_x_midpoint_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_x_midpoint_ptr.push_back(&F_x_midpoint[ei][idx_face_x]);
                                }
                                
                                d_Riemann_solver.computeIntercellFluxForSingleSpecies(
                                    F_x_midpoint_ptr,
                                    &rho_L,
                                    &rho_R,
                                    m_L_ptr,
                                    m_R_ptr,
                                    &E_L,
                                    &E_R,
                                    X_DIRECTION);
                            }
                        }
                    }
                    
                    // Compute the mid-point flux in the y direction.
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = -1; j < interior_dims[1] + 2; j++)
                            {
                                // Compute the index of face of mid-point flux and
                                // projection matrix.
                                const int idx_face_y = (j + 1) +
                                    (k + 1)*(interior_dims[1] + 3) +
                                    (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                
                                boost::multi_array<double, 2> W_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                /*
                                 * Project primitive variables onto characteristic fields.
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
                                    
                                    std::vector<const double*> V_ptr;
                                    V_ptr.push_back(&rho[idx_cell]);
                                    V_ptr.push_back(&u[idx_cell]);
                                    V_ptr.push_back(&v[idx_cell]);
                                    V_ptr.push_back(&w[idx_cell]);
                                    V_ptr.push_back(&p[idx_cell]);
                                    
                                    std::vector<double*> W_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        W_ptr.push_back(&W_array[m][ei]);
                                    }
                                    
                                    projectPrimitiveVariablesToCharacteristicFields(
                                        W_ptr,
                                        V_ptr,
                                        R_y_inv_intercell);
                                }
                                
                                /*
                                 * Do WENO interplation on characteristic variables to get W_B
                                 * and W_T
                                 */
                                
                                std::vector<double> W_B;
                                std::vector<double> W_T;
                                
                                performWENOInterpolation(W_B, W_T, W_array);
                                
                                /*
                                 * Project characteristic variables back to physical fields.
                                 */
                                
                                double rho_B;
                                double rho_T;
                                
                                std::vector<double> vel_B;
                                std::vector<double> vel_T;
                                vel_B.resize(d_dim.getValue());
                                vel_T.resize(d_dim.getValue());
                                
                                double p_B;
                                double p_T;
                                
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
                                
                                std::vector<const double*> W_B_ptr;
                                std::vector<const double*> W_T_ptr;
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_B_ptr.push_back(&W_B[ei]);
                                    W_T_ptr.push_back(&W_T[ei]);
                                }
                                
                                std::vector<double*> V_B_ptr;
                                V_B_ptr.push_back(&rho_B);
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_B_ptr.push_back(&vel_B[di]);
                                }
                                V_B_ptr.push_back(&p_B);
                                
                                std::vector<double*> V_T_ptr;
                                V_T_ptr.push_back(&rho_T);
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_T_ptr.push_back(&vel_T[di]);
                                }
                                V_T_ptr.push_back(&p_T);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_B_ptr,
                                    W_B_ptr,
                                    R_y_intercell);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_T_ptr,
                                    W_T_ptr,
                                    R_y_intercell);
                                
                                /*
                                 * Convert the primitive variables into conservative variables.
                                 */
                                
                                std::vector<double> m_B;
                                std::vector<double> m_T;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_B.push_back(rho_B*vel_B[di]);
                                    m_T.push_back(rho_T*vel_T[di]);
                                }
                                
                                std::vector<const double*> vel_B_ptr;
                                std::vector<const double*> vel_T_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    vel_B_ptr.push_back(&vel_B[di]);
                                    vel_T_ptr.push_back(&vel_T[di]);
                                }
                                
                                double E_B = d_equation_of_state->
                                    getTotalEnergy(
                                        &rho_B,
                                        vel_B_ptr,
                                        &p_B);
                                
                                double E_T = d_equation_of_state->
                                    getTotalEnergy(
                                        &rho_T,
                                        vel_T_ptr,
                                        &p_T);
                                
                                /*
                                 * If the WENO interpolated density or energy are negative,
                                 * use constant interpolation.
                                 */
                                if ((rho_B < 0) || (rho_T < 0) || (E_B < 0) || (E_T < 0))
                                {
                                    // Compute the indices of bottom cell and top cell.
                                    const int idx_cell_B = (i + d_num_ghosts[0]) +
                                        (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_T = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    rho_B = rho[idx_cell_B];
                                    rho_T = rho[idx_cell_T];
                                    
                                    m_B[0] = rho_u[idx_cell_B];
                                    m_B[1] = rho_v[idx_cell_B];
                                    m_B[2] = rho_w[idx_cell_B];
                                    m_T[0] = rho_u[idx_cell_T];
                                    m_T[1] = rho_v[idx_cell_T];
                                    m_T[2] = rho_w[idx_cell_T];
                                    
                                    E_B = E[idx_cell_B];
                                    E_T = E[idx_cell_T];
                                }
                                
                                /*
                                 * Apply the Riemann solver.
                                 */
                                
                                std::vector<const double*> m_B_ptr;
                                std::vector<const double*> m_T_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_B_ptr.push_back(&m_B[di]);
                                    m_T_ptr.push_back(&m_T[di]);
                                }
                                
                                std::vector<double*> F_y_midpoint_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_y_midpoint_ptr.push_back(&F_y_midpoint[ei][idx_face_y]);
                                }
                                
                                d_Riemann_solver.computeIntercellFluxForSingleSpecies(
                                    F_y_midpoint_ptr,
                                    &rho_B,
                                    &rho_T,
                                    m_B_ptr,
                                    m_T_ptr,
                                    &E_B,
                                    &E_T,
                                    Y_DIRECTION);
                            }
                        }
                    }
                    
                    // Compute the mid-point flux in the z direction.
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int k = -1; k < interior_dims[2] + 2; k++)
                            {
                                // Compute the index of face of mid-point flux and
                                // projection matrix.
                                const int idx_face_z = (k + 1) +
                                    (i + 1)*(interior_dims[2] + 3) +
                                    (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                
                                boost::multi_array<double, 2> W_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                /*
                                 * Project primitive variables onto characteristic fields.
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
                                    
                                    std::vector<const double*> V_ptr;
                                    V_ptr.push_back(&rho[idx_cell]);
                                    V_ptr.push_back(&u[idx_cell]);
                                    V_ptr.push_back(&v[idx_cell]);
                                    V_ptr.push_back(&w[idx_cell]);
                                    V_ptr.push_back(&p[idx_cell]);
                                    
                                    std::vector<double*> W_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        W_ptr.push_back(&W_array[m][ei]);
                                    }
                                    
                                    projectPrimitiveVariablesToCharacteristicFields(
                                        W_ptr,
                                        V_ptr,
                                        R_z_inv_intercell);
                                }
                                
                                /*
                                 * Do WENO interplation on characteristic variables to get W_B
                                 * and W_F
                                 */
                                
                                std::vector<double> W_B;
                                std::vector<double> W_F;
                                
                                performWENOInterpolation(W_B, W_F, W_array);
                                
                                /*
                                 * Project characteristic variables back to physical fields.
                                 */
                                
                                double rho_B;
                                double rho_F;
                                
                                std::vector<double> vel_B;
                                std::vector<double> vel_F;
                                vel_B.resize(d_dim.getValue());
                                vel_F.resize(d_dim.getValue());
                                
                                double p_B;
                                double p_F;
                                
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
                                
                                std::vector<const double*> W_B_ptr;
                                std::vector<const double*> W_F_ptr;
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_B_ptr.push_back(&W_B[ei]);
                                    W_F_ptr.push_back(&W_F[ei]);
                                }
                                
                                std::vector<double*> V_B_ptr;
                                V_B_ptr.push_back(&rho_B);
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_B_ptr.push_back(&vel_B[di]);
                                }
                                V_B_ptr.push_back(&p_B);
                                
                                std::vector<double*> V_F_ptr;
                                V_F_ptr.push_back(&rho_F);
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_F_ptr.push_back(&vel_F[di]);
                                }
                                V_F_ptr.push_back(&p_F);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_B_ptr,
                                    W_B_ptr,
                                    R_z_intercell);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_F_ptr,
                                    W_F_ptr,
                                    R_z_intercell);
                                
                                /*
                                 * Convert the primitive variables into conservative variables.
                                 */
                                
                                std::vector<double> m_B;
                                std::vector<double> m_F;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_B.push_back(rho_B*vel_B[di]);
                                    m_F.push_back(rho_F*vel_F[di]);
                                }
                                
                                std::vector<const double*> vel_B_ptr;
                                std::vector<const double*> vel_F_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    vel_B_ptr.push_back(&vel_B[di]);
                                    vel_F_ptr.push_back(&vel_F[di]);
                                }
                                
                                double E_B = d_equation_of_state->
                                    getTotalEnergy(
                                        &rho_B,
                                        vel_B_ptr,
                                        &p_B);
                                
                                double E_F = d_equation_of_state->
                                    getTotalEnergy(
                                        &rho_F,
                                        vel_F_ptr,
                                        &p_F);
                                
                                /*
                                 * If the WENO interpolated density or energy are negative,
                                 * use constant interpolation.
                                 */
                                if ((rho_B < 0) || (rho_F < 0) || (E_B < 0) || (E_F < 0))
                                {
                                    // Compute the indices of back cell and front cell.
                                    const int idx_cell_B = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_F = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                        
                                    rho_B = rho[idx_cell_B];
                                    rho_F = rho[idx_cell_F];
                                    
                                    m_B[0] = rho_u[idx_cell_B];
                                    m_B[1] = rho_v[idx_cell_B];
                                    m_B[2] = rho_w[idx_cell_B];
                                    m_F[0] = rho_u[idx_cell_F];
                                    m_F[1] = rho_v[idx_cell_F];
                                    m_F[2] = rho_w[idx_cell_F];
                                    
                                    E_B = E[idx_cell_B];
                                    E_F = E[idx_cell_F];
                                }
                                
                                /*
                                 * Apply the Riemann solver.
                                 */
                                
                                std::vector<const double*> m_B_ptr;
                                std::vector<const double*> m_F_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_B_ptr.push_back(&m_B[di]);
                                    m_F_ptr.push_back(&m_F[di]);
                                }
                                
                                std::vector<double*> F_z_midpoint_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_z_midpoint_ptr.push_back(&F_z_midpoint[ei][idx_face_z]);
                                }
                                
                                d_Riemann_solver.computeIntercellFluxForSingleSpecies(
                                    F_z_midpoint_ptr,
                                    &rho_B,
                                    &rho_F,
                                    m_B_ptr,
                                    m_F_ptr,
                                    &E_B,
                                    &E_F,
                                    Z_DIRECTION);
                            }
                        }
                    }
                    
                    // Compute the flux in the x direction.
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0] + 1; i++)
                            {
                                // Compute the indices.
                                const int idx_face_x = i +
                                    j*(interior_dims[0] + 1) +
                                    k*(interior_dims[0] + 1)*interior_dims[1];
                                
                                const int idx_midpoint_x = (i + 1) +
                                    (j + 1)*(interior_dims[0] + 3) +
                                    (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                
                                const int idx_node_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_node_R = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                // Compute the flux.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    convective_flux->getPointer(0, ei)[idx_face_x] = dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x + 1] +
                                        F_x_midpoint[ei][idx_midpoint_x - 1]) -
                                        3.0/10*(F_x_node[ei][idx_node_R] +
                                        F_x_node[ei][idx_node_L]) +
                                        23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
                                }
                            }
                        }
                    }
                    
                    // Compute the flux in the y direction.
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = 0; j < interior_dims[1] + 1; j++)
                            {
                                // Compute the indices.
                                const int idx_face_y = j +
                                    k*(interior_dims[1] + 1) +
                                    i*(interior_dims[1] + 1)*interior_dims[2];
                                
                                const int idx_midpoint_y = (j + 1) +
                                    (k + 1)*(interior_dims[1] + 3) +
                                    (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                
                                const int idx_node_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_node_T = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                // Compute the flux.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    convective_flux->getPointer(1, ei)[idx_face_y] =
                                        dt*(1.0/30*(F_y_midpoint[ei][idx_midpoint_y + 1] +
                                        F_y_midpoint[ei][idx_midpoint_y - 1]) -
                                        3.0/10*(F_y_node[ei][idx_node_T] +
                                        F_y_node[ei][idx_node_B]) +
                                        23.0/15*F_y_midpoint[ei][idx_midpoint_y]);
                                }
                            }
                        }
                    }
                    
                    // Compute the flux in the z direction.
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int k = 0; k < interior_dims[2] + 1; k++)
                            {
                                // Compute the indices.
                                const int idx_face_z = k +
                                    i*(interior_dims[2] + 1) +
                                    j*(interior_dims[2] + 1)*interior_dims[0];
                                
                                const int idx_midpoint_z = (k + 1) +
                                    (i + 1)*(interior_dims[2] + 3) +
                                    (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                
                                const int idx_node_B = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_node_F = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                // Compute the flux.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    convective_flux->getPointer(2, ei)[idx_face_z] =
                                        dt*(1.0/30*(F_z_midpoint[ei][idx_midpoint_z + 1] +
                                        F_z_midpoint[ei][idx_midpoint_z - 1]) -
                                        3.0/10*(F_z_node[ei][idx_node_F] +
                                        F_z_node[ei][idx_node_B]) +
                                        23.0/15*F_z_midpoint[ei][idx_midpoint_z]);
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
                    new pdat::FaceData<double>(interior_box, d_dim.getValue(), hier::IntVector::getOne(d_dim)));
                
                boost::shared_ptr<pdat::CellData<double> > velocity(
                    new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
                
                boost::shared_ptr<pdat::CellData<double> > pressure(
                    new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                
                boost::shared_ptr<pdat::CellData<double> > sound_speed(
                    new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > convective_flux_node;
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    convective_flux_node.push_back(boost::make_shared<pdat::CellData<double> >(
                        interior_box, d_num_eqn, d_num_ghosts));
                }
                
                boost::shared_ptr<pdat::FaceData<double> > convective_flux_midpoint(
                new pdat::FaceData<double>(interior_box, d_num_eqn, hier::IntVector::getOne(d_dim)));
                
                boost::shared_ptr<pdat::FaceData<double> > projection_matrix(
                    new pdat::FaceData<double>(
                        interior_box,
                        d_num_eqn*d_num_eqn,
                        hier::IntVector::getOne(d_dim)));
                
                boost::shared_ptr<pdat::FaceData<double> > projection_matrix_inv(
                    new pdat::FaceData<double>(
                        interior_box,
                        d_num_eqn*d_num_eqn,
                        hier::IntVector::getOne(d_dim)));
                
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
                    std::vector<double*> F_x_node;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
                    }
                    std::vector<double*> F_x_midpoint;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_midpoint.push_back(convective_flux_midpoint->getPointer(0, ei));
                    }
                    
                    // Compute the field of velocities, pressure and fluxes.
                    for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                    {
                        // Compute index into linear data array.
                        const int idx = i + d_num_ghosts[0];
                        
                        std::vector<const double*> Y_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Y_ptr.push_back(&Y[si][idx]);
                        }
                        
                        u[idx] = rho_u[idx]/rho[idx];
                        
                        std::vector<const double*> m_ptr;
                        m_ptr.push_back(&rho_u[idx]);
                        
                        p[idx] = d_equation_of_state->
                            getPressureWithMassFraction(
                                &rho[idx],
                                m_ptr,
                                &E[idx],
                                Y_ptr);
                        
                        c[idx] = d_equation_of_state->
                            getSoundSpeedWithMassFractionAndPressure(
                                &rho[idx],
                                Y_ptr,
                                &p[idx]);
                        
                        F_x_node[0][idx] = rho_u[idx];
                        F_x_node[1][idx] = rho_u[idx]*u[idx] + p[idx];
                        F_x_node[2][idx] = u[idx]*(E[idx] + p[idx]);
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            F_x_node[3 + si][idx] = u[idx]*Y[si][idx];
                        }
                    }
                    
                    /*
                     * Compute the projection matrix and its inverse at the face normal to the
                     * x direction.
                     */
                    for (int i = -1; i < interior_dims[0] + 2; i++)
                    {
                        // Compute the indices of left cell, right cell and face of
                        // projection matrix.
                        const int idx_cell_L = i - 1 + d_num_ghosts[0];
                        
                        const int idx_cell_R = i + d_num_ghosts[0];
                        
                        const int idx_face_x = i + 1;
                        
                        // Get left and right quantities.
                        const double& rho_L = rho[idx_cell_L];
                        const double& rho_R = rho[idx_cell_R];
                        
                        const double& c_L = c[idx_cell_L];
                        const double& c_R = c[idx_cell_R];
                        
                        // Compute simply-averaged quantities.
                        const double rho_average = 0.5*(rho_L + rho_R);
                        const double c_average = 0.5*(c_L + c_R);
                        
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
                                    
                                *R_x_intercell[ei][ej] = 0.0;
                                *R_x_inv_intercell[ei][ej] = 0.0;
                            }
                        }
                        
                        /*
                         * Compute R_x_intercell.
                         */
                        
                        *R_x_intercell[0][0] = 1.0/(c_average*c_average);
                        *R_x_intercell[0][1] = 1.0;
                        *R_x_intercell[0][d_num_eqn - 1] = 1.0/(c_average*c_average);
                        
                        *R_x_intercell[1][0] = -1.0/(rho_average*c_average);
                        *R_x_intercell[1][d_num_eqn - 1] = 1.0/(rho_average*c_average);
                        
                        *R_x_intercell[2][0] = 1.0;
                        *R_x_intercell[2][d_num_eqn - 1] = 1.0;
                        
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            *R_x_intercell[3 + si][2 + si] = 1.0;
                        }
                        
                        /*
                         * Compute R_x_inv_intercell.
                         */
                        
                        *R_x_inv_intercell[0][1] = -0.5*rho_average*c_average;
                        *R_x_inv_intercell[0][2] = 0.5;
                        
                        *R_x_inv_intercell[1][0] = 1.0;
                        *R_x_inv_intercell[1][2] = -1.0/(c_average*c_average);
                        
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            *R_x_inv_intercell[2 + si][3 + si] = 1.0;
                        }
                        
                        *R_x_inv_intercell[d_num_species + 1][1] = 0.5*rho_average*c_average;
                        *R_x_inv_intercell[d_num_species + 1][2] = 0.5;
                    }
                    
                    // Compute the mid-point flux in the x direction.
                    for (int i = -1; i < interior_dims[0] + 2; i++)
                    {
                        // Compute the index of face of mid-point flux and
                        // projection matrix.
                        const int idx_face_x = i + 1;
                        
                        boost::multi_array<double, 2> W_array(
                            boost::extents[6][d_num_eqn]);
                        
                        /*
                         * Project primitive variables onto characteristic fields.
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
                            
                            std::vector<const double*> V_ptr;
                            V_ptr.push_back(&rho[idx_cell]);
                            V_ptr.push_back(&u[idx_cell]);
                            V_ptr.push_back(&p[idx_cell]);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                V_ptr.push_back(&Y[si][idx_cell]);
                            }
                            
                            std::vector<double*> W_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                W_ptr.push_back(&W_array[m][ei]);
                            }
                            projectPrimitiveVariablesToCharacteristicFields(
                                W_ptr,
                                V_ptr,
                                R_x_inv_intercell);
                        }
                        
                        /*
                         * Do WENO interplation on characteristic variables to get W_L
                         * and W_R
                         */
                        
                        std::vector<double> W_L;
                        std::vector<double> W_R;
                        
                        performWENOInterpolation(W_L, W_R, W_array);
                        
                        /*
                         * Project characteristic variables back to physical fields.
                         */
                        
                        double rho_L;
                        double rho_R;
                        
                        std::vector<double> vel_L;
                        std::vector<double> vel_R;
                        vel_L.resize(d_dim.getValue());
                        vel_R.resize(d_dim.getValue());
                        
                        double p_L;
                        double p_R;
                        
                        std::vector<double> Y_L;
                        std::vector<double> Y_R;
                        Y_L.resize(d_num_species - 1);
                        Y_R.resize(d_num_species - 1);
                        
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
                        
                        std::vector<const double*> W_L_ptr;
                        std::vector<const double*> W_R_ptr;
                        
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            W_L_ptr.push_back(&W_L[ei]);
                            W_R_ptr.push_back(&W_R[ei]);
                        }
                        
                        std::vector<double*> V_L_ptr;
                        V_L_ptr.push_back(&rho_L);
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            V_L_ptr.push_back(&vel_L[di]);
                        }
                        V_L_ptr.push_back(&p_L);
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            V_L_ptr.push_back(&Y_L[si]);
                        }
                        
                        std::vector<double*> V_R_ptr;
                        V_R_ptr.push_back(&rho_R);
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            V_R_ptr.push_back(&vel_R[di]);
                        }
                        V_R_ptr.push_back(&p_R);
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            V_R_ptr.push_back(&Y_R[si]);
                        }
                        
                        projectCharacteristicVariablesToPhysicalFields(
                            V_L_ptr,
                            W_L_ptr,
                            R_x_intercell);
                        
                        projectCharacteristicVariablesToPhysicalFields(
                            V_R_ptr,
                            W_R_ptr,
                            R_x_intercell);
                        
                        /*
                         * Convert the primitive variables into conservative variables.
                         */
                        
                        std::vector<double> m_L;
                        std::vector<double> m_R;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            m_L.push_back(rho_L*vel_L[di]);
                            m_R.push_back(rho_R*vel_R[di]);
                        }
                        
                        std::vector<const double*> vel_L_ptr;
                        std::vector<const double*> vel_R_ptr;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            vel_L_ptr.push_back(&vel_L[di]);
                            vel_R_ptr.push_back(&vel_R[di]);
                        }
                        
                        std::vector<const double*> Y_L_ptr;
                        std::vector<const double*> Y_R_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Y_L_ptr.push_back(&Y_L[si]);
                            Y_R_ptr.push_back(&Y_R[si]);
                        }
                        
                        double E_L = d_equation_of_state->
                            getTotalEnergyWithMassFraction(
                                &rho_L,
                                vel_L_ptr,
                                &p_L,
                                Y_L_ptr);
                        
                        double E_R = d_equation_of_state->
                            getTotalEnergyWithMassFraction(
                                &rho_R,
                                vel_R_ptr,
                                &p_R,
                                Y_R_ptr);
                        
                        /*
                         * If the WENO interpolated density or energy are negative,
                         * use constant interpolation.
                         */
                        if ((rho_L < 0) || (rho_R < 0) || (E_L < 0) || (E_R < 0))
                        {
                            // Compute the indices of left cell and right cell.
                            const int idx_cell_L = i - 1 + d_num_ghosts[0];
                            const int idx_cell_R = i + d_num_ghosts[0];
                            
                            rho_L = rho[idx_cell_L];
                            rho_R = rho[idx_cell_R];
                            
                            m_L[0] = rho_u[idx_cell_L];
                            m_R[0] = rho_u[idx_cell_R];
                            
                            E_L = E[idx_cell_L];
                            E_R = E[idx_cell_R];
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Y_L[si] = Y[si][idx_cell_L];
                                Y_R[si] = Y[si][idx_cell_R];
                            }
                        }
                        
                        /*
                         * Apply the Riemann solver.
                         */
                        
                        std::vector<const double*> m_L_ptr;
                        std::vector<const double*> m_R_ptr;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            m_L_ptr.push_back(&m_L[di]);
                            m_R_ptr.push_back(&m_R[di]);
                        }
                        
                        std::vector<double*> F_x_midpoint_ptr;
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            F_x_midpoint_ptr.push_back(&F_x_midpoint[ei][idx_face_x]);
                        }
                        
                        std::vector<double*> vel_x_ptr;
                        for (int vi = 0; vi < d_dim.getValue(); vi++)
                        {
                            vel_x_ptr.push_back(&(velocity_intercell->getPointer(0, vi)[idx_face_x]));
                        }
                        
                        d_Riemann_solver.computeIntercellFluxAndVelocityForFourEqnShyue(
                            F_x_midpoint_ptr,
                            vel_x_ptr,
                            &rho_L,
                            &rho_R,
                            m_L_ptr,
                            m_R_ptr,
                            &E_L,
                            &E_R,
                            Y_L_ptr,
                            Y_R_ptr,
                            X_DIRECTION);
                    }
                    
                    // Compute the flux in the x direction.
                    for (int i = 0; i < interior_dims[0] + 1; i++)
                    {
                        // Compute the indices.
                        const int idx_face_x = i;
                        
                        const int idx_midpoint_x = i + 1;
                        
                        const int idx_node_L = i - 1 + d_num_ghosts[0];
                        
                        const int idx_node_R = i + d_num_ghosts[0];
                        
                        // Compute the flux.
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            convective_flux->getPointer(0, ei)[idx_face_x] = dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x + 1] +
                                F_x_midpoint[ei][idx_midpoint_x - 1]) -
                                3.0/10*(F_x_node[ei][idx_node_R] +
                                F_x_node[ei][idx_node_L]) +
                                23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
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
                            
                            const int idx_face_x_L = i + 1;
                            
                            const int idx_face_x_R = i + 2;
                            
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
                    std::vector<double*> F_x_node;
                    std::vector<double*> F_y_node;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
                        F_y_node.push_back(convective_flux_node[1]->getPointer(ei));
                    }
                    std::vector<double*> F_x_midpoint;
                    std::vector<double*> F_y_midpoint;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_midpoint.push_back(convective_flux_midpoint->getPointer(0, ei));
                        F_y_midpoint.push_back(convective_flux_midpoint->getPointer(1, ei));
                    }
                    
                    // Compute the field of velocities, pressure and fluxes.
                    for (int j = -d_num_ghosts[1]; j < interior_dims[1] + d_num_ghosts[1]; j++)
                    {
                        for (int i = -d_num_ghosts[0]; i < interior_dims[0] + d_num_ghosts[0]; i++)
                        {
                            // Compute index into linear data array.
                            const int idx = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            std::vector<const double*> Y_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Y_ptr.push_back(&Y[si][idx]);
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
                                    Y_ptr);
                            
                            c[idx] = d_equation_of_state->
                                getSoundSpeedWithMassFractionAndPressure(
                                    &rho[idx],
                                    Y_ptr,
                                    &p[idx]);
                            
                            F_x_node[0][idx] = rho_u[idx];
                            F_x_node[1][idx] = rho_u[idx]*u[idx] + p[idx];
                            F_x_node[2][idx] = rho_u[idx]*v[idx];
                            F_x_node[3][idx] = u[idx]*(E[idx] + p[idx]);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_x_node[4 + si][idx] = u[idx]*Y[si][idx];
                            }
                            
                            F_y_node[0][idx] = rho_v[idx];
                            F_y_node[1][idx] = rho_v[idx]*u[idx];
                            F_y_node[2][idx] = rho_v[idx]*v[idx] + p[idx];
                            F_y_node[3][idx] = v[idx]*(E[idx] + p[idx]);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_y_node[4 + si][idx] = v[idx]*Y[si][idx];
                            }
                        }
                    }
                    
                    /*
                     * Compute the projection matrix and its inverse at the face normal to the
                     * x direction.
                     */
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = -1; i < interior_dims[0] + 2; i++)
                        {
                            // Compute the indices of left cell, right cell and face of
                            // projection matrix.
                            const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_cell_R = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_face_x = (i + 1) +
                                (j + 1)*(interior_dims[0] + 3);
                            
                            // Get left and right quantities.
                            const double& rho_L = rho[idx_cell_L];
                            const double& rho_R = rho[idx_cell_R];
                            
                            const double& c_L = c[idx_cell_L];
                            const double& c_R = c[idx_cell_R];
                            
                            // Compute simply-averaged quantities.
                            const double rho_average = 0.5*(rho_L + rho_R);
                            const double c_average = 0.5*(c_L + c_R);
                            
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
                                        
                                    *R_x_intercell[ei][ej] = 0.0;
                                    *R_x_inv_intercell[ei][ej] = 0.0;
                                }
                            }
                            
                            /*
                             * Compute R_x_intercell.
                             */
                            
                            *R_x_intercell[0][0] = 1.0/(c_average*c_average);
                            *R_x_intercell[0][1] = 1.0;
                            *R_x_intercell[0][d_num_eqn - 1] = 1.0/(c_average*c_average);
                            
                            *R_x_intercell[1][0] = -1.0/(rho_average*c_average);
                            *R_x_intercell[1][d_num_eqn - 1] = 1.0/(rho_average*c_average);
                            
                            *R_x_intercell[2][2] = 1.0;
                            
                            *R_x_intercell[3][0] = 1.0;
                            *R_x_intercell[3][d_num_eqn - 1] = 1.0;
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                *R_x_intercell[4 + si][3 + si] = 1.0;
                            }
                            
                            /*
                             * Compute R_x_inv_intercell.
                             */
                            
                            *R_x_inv_intercell[0][1] = -0.5*rho_average*c_average;
                            *R_x_inv_intercell[0][3] = 0.5;
                            
                            *R_x_inv_intercell[1][0] = 1.0;
                            *R_x_inv_intercell[1][3] = -1.0/(c_average*c_average);
                            
                            *R_x_inv_intercell[2][2] = 1.0;
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                *R_x_inv_intercell[3 + si][4 + si] = 1.0;
                            }
                            
                            *R_x_inv_intercell[d_num_species + 2][1] = 0.5*rho_average*c_average;
                            *R_x_inv_intercell[d_num_species + 2][3] = 0.5;
                        }
                    }
                    
                    /*
                     * Compute the projection matrix and its inverse at the face normal to the
                     * y direction.
                     */
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int j = -1; j < interior_dims[1] + 2; j++)
                        {
                            // Compute the indices of bottom cell, top cell and face of
                            // projection matrix.
                            const int idx_cell_B = (i + d_num_ghosts[0]) +
                                (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_cell_T = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_face_y = (j + 1) +
                                (i + 1)*(interior_dims[1] + 3);
                            
                            // Get bottom and top quantities.
                            
                            const double& rho_B = rho[idx_cell_B];
                            const double& rho_T = rho[idx_cell_T];
                            
                            const double& c_B = c[idx_cell_B];
                            const double& c_T = c[idx_cell_T];
                            
                            // Compute simply-averaged quantities.
                            const double rho_average = 0.5*(rho_B + rho_T);
                            const double c_average = 0.5*(c_B + c_T);
                            
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
                                    
                                    *R_y_intercell[ei][ej] = 0.0;
                                    *R_y_inv_intercell[ei][ej] = 0.0;
                                }
                            }
                            
                            /*
                             * Compute R_y_intercell.
                             */
                            
                            *R_y_intercell[0][0] = 1.0/(c_average*c_average);
                            *R_y_intercell[0][1] = 1.0;
                            *R_y_intercell[0][d_num_eqn - 1] = 1.0/(c_average*c_average);
                            
                            *R_y_intercell[1][2] = 1.0;
                            
                            *R_y_intercell[2][0] = -1.0/(rho_average*c_average);
                            *R_y_intercell[2][d_num_eqn - 1] = 1.0/(rho_average*c_average);
                            
                            *R_y_intercell[3][0] = 1.0;
                            *R_y_intercell[3][d_num_eqn - 1] = 1.0;
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                *R_y_intercell[4 + si][3 + si] = 1.0;
                            }
                            
                            /*
                             * Compute R_y_inv_intercell.
                             */
                            
                            *R_y_inv_intercell[0][2] = -0.5*rho_average*c_average;
                            *R_y_inv_intercell[0][3] = 0.5;
                            
                            *R_y_inv_intercell[1][0] = 1.0;
                            *R_y_inv_intercell[1][3] = -1.0/(c_average*c_average);
                            
                            *R_y_inv_intercell[2][1] = 1.0;
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                *R_y_inv_intercell[3 + si][4 + si] = 1.0;
                            }
                            
                            *R_y_inv_intercell[d_num_species + 2][2] = 0.5*rho_average*c_average;
                            *R_y_inv_intercell[d_num_species + 2][3] = 0.5;
                        }
                    }
                    
                    // Compute the mid-point flux in the x direction.
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = -1; i < interior_dims[0] + 2; i++)
                        {
                            // Compute the index of face of mid-point flux and
                            // projection matrix.
                            const int idx_face_x = (i + 1) +
                                (j + 1)*(interior_dims[0] + 3);
                            
                            boost::multi_array<double, 2> W_array(
                                boost::extents[6][d_num_eqn]);
                            
                            /*
                             * Project primitive variables onto characteristic fields.
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
                                
                                std::vector<const double*> V_ptr;
                                V_ptr.push_back(&rho[idx_cell]);
                                V_ptr.push_back(&u[idx_cell]);
                                V_ptr.push_back(&v[idx_cell]);
                                V_ptr.push_back(&p[idx_cell]);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_ptr.push_back(&Y[si][idx_cell]);
                                }
                                
                                std::vector<double*> W_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_ptr.push_back(&W_array[m][ei]);
                                }
                                projectPrimitiveVariablesToCharacteristicFields(
                                    W_ptr,
                                    V_ptr,
                                    R_x_inv_intercell);
                            }
                            
                            /*
                             * Do WENO interplation on characteristic variables to get W_L
                             * and W_R
                             */
                            
                            std::vector<double> W_L;
                            std::vector<double> W_R;
                            
                            performWENOInterpolation(W_L, W_R, W_array);
                            
                            /*
                             * Project characteristic variables back to physical fields.
                             */
                            
                            double rho_L;
                            double rho_R;
                            
                            std::vector<double> vel_L;
                            std::vector<double> vel_R;
                            vel_L.resize(d_dim.getValue());
                            vel_R.resize(d_dim.getValue());
                            
                            double p_L;
                            double p_R;
                            
                            std::vector<double> Y_L;
                            std::vector<double> Y_R;
                            Y_L.resize(d_num_species - 1);
                            Y_R.resize(d_num_species - 1);
                            
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
                            
                            std::vector<const double*> W_L_ptr;
                            std::vector<const double*> W_R_ptr;
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                W_L_ptr.push_back(&W_L[ei]);
                                W_R_ptr.push_back(&W_R[ei]);
                            }
                            
                            std::vector<double*> V_L_ptr;
                            V_L_ptr.push_back(&rho_L);
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                V_L_ptr.push_back(&vel_L[di]);
                            }
                            V_L_ptr.push_back(&p_L);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                V_L_ptr.push_back(&Y_L[si]);
                            }
                            
                            std::vector<double*> V_R_ptr;
                            V_R_ptr.push_back(&rho_R);
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                V_R_ptr.push_back(&vel_R[di]);
                            }
                            V_R_ptr.push_back(&p_R);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                V_R_ptr.push_back(&Y_R[si]);
                            }
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                V_L_ptr,
                                W_L_ptr,
                                R_x_intercell);
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                V_R_ptr,
                                W_R_ptr,
                                R_x_intercell);
                            
                            /*
                             * Convert the primitive variables into conservative variables.
                             */
                            
                            std::vector<double> m_L;
                            std::vector<double> m_R;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_L.push_back(rho_L*vel_L[di]);
                                m_R.push_back(rho_R*vel_R[di]);
                            }
                            
                            std::vector<const double*> vel_L_ptr;
                            std::vector<const double*> vel_R_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                vel_L_ptr.push_back(&vel_L[di]);
                                vel_R_ptr.push_back(&vel_R[di]);
                            }
                            
                            std::vector<const double*> Y_L_ptr;
                            std::vector<const double*> Y_R_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Y_L_ptr.push_back(&Y_L[si]);
                                Y_R_ptr.push_back(&Y_R[si]);
                            }
                            
                            double E_L = d_equation_of_state->
                                getTotalEnergyWithMassFraction(
                                    &rho_L,
                                    vel_L_ptr,
                                    &p_L,
                                    Y_L_ptr);
                            
                            double E_R = d_equation_of_state->
                                getTotalEnergyWithMassFraction(
                                    &rho_R,
                                    vel_R_ptr,
                                    &p_R,
                                    Y_R_ptr);
                            
                            /*
                             * If the WENO interpolated density or energy are negative,
                             * use constant interpolation.
                             */
                            if ((rho_L < 0) || (rho_R < 0) || (E_L < 0) || (E_R < 0))
                            {
                                // Compute the indices of left cell and right cell.
                                const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_R = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                rho_L = rho[idx_cell_L];
                                rho_R = rho[idx_cell_R];
                                
                                m_L[0] = rho_u[idx_cell_L];
                                m_L[1] = rho_v[idx_cell_L];
                                m_R[0] = rho_u[idx_cell_R];
                                m_R[1] = rho_v[idx_cell_R];
                                
                                E_L = E[idx_cell_L];
                                E_R = E[idx_cell_R];
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_L[si] = Y[si][idx_cell_L];
                                    Y_R[si] = Y[si][idx_cell_R];
                                }
                            }
                            
                            /*
                             * Apply the Riemann solver.
                             */
                            
                            std::vector<const double*> m_L_ptr;
                            std::vector<const double*> m_R_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_L_ptr.push_back(&m_L[di]);
                                m_R_ptr.push_back(&m_R[di]);
                            }
                            
                            std::vector<double*> F_x_midpoint_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_x_midpoint_ptr.push_back(&F_x_midpoint[ei][idx_face_x]);
                            }
                            
                            std::vector<double*> vel_x_ptr;
                            for (int vi = 0; vi < d_dim.getValue(); vi++)
                            {
                                vel_x_ptr.push_back(&(velocity_intercell->getPointer(0, vi)[idx_face_x]));
                            }
                            
                            d_Riemann_solver.computeIntercellFluxAndVelocityForFourEqnShyue(
                                F_x_midpoint_ptr,
                                vel_x_ptr,
                                &rho_L,
                                &rho_R,
                                m_L_ptr,
                                m_R_ptr,
                                &E_L,
                                &E_R,
                                Y_L_ptr,
                                Y_R_ptr,
                                X_DIRECTION);
                        }
                    }
                    
                    // Compute the mid-point flux in the y direction.
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int j = -1; j < interior_dims[1] + 2; j++)
                        {
                            // Compute the index of face of mid-point flux and
                            // projection matrix.
                            const int idx_face_y = (j + 1) +
                                (i + 1)*(interior_dims[1] + 3);
                            
                            boost::multi_array<double, 2> W_array(
                                boost::extents[6][d_num_eqn]);
                            
                            /*
                             * Project primitive variables onto characteristic fields.
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
                                
                                std::vector<const double*> V_ptr;
                                V_ptr.push_back(&rho[idx_cell]);
                                V_ptr.push_back(&u[idx_cell]);
                                V_ptr.push_back(&v[idx_cell]);
                                V_ptr.push_back(&p[idx_cell]);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_ptr.push_back(&Y[si][idx_cell]);
                                }
                                
                                std::vector<double*> W_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_ptr.push_back(&W_array[m][ei]);
                                }
                                
                                projectPrimitiveVariablesToCharacteristicFields(
                                    W_ptr,
                                    V_ptr,
                                    R_y_inv_intercell);
                            }
                            
                            /*
                             * Do WENO interplation on characteristic variables to get W_B
                             * and W_T
                             */
                            
                            std::vector<double> W_B;
                            std::vector<double> W_T;
                            
                            performWENOInterpolation(W_B, W_T, W_array);
                            
                            /*
                             * Project characteristic variables back to physical fields.
                             */
                            
                            double rho_B;
                            double rho_T;
                            
                            std::vector<double> vel_B;
                            std::vector<double> vel_T;
                            vel_B.resize(d_dim.getValue());
                            vel_T.resize(d_dim.getValue());
                            
                            double p_B;
                            double p_T;
                            
                            std::vector<double> Y_B;
                            std::vector<double> Y_T;
                            Y_B.resize(d_num_species - 1);
                            Y_T.resize(d_num_species - 1);
                            
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
                            
                            std::vector<const double*> W_B_ptr;
                            std::vector<const double*> W_T_ptr;
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                W_B_ptr.push_back(&W_B[ei]);
                                W_T_ptr.push_back(&W_T[ei]);
                            }
                            
                            std::vector<double*> V_B_ptr;
                            V_B_ptr.push_back(&rho_B);
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                V_B_ptr.push_back(&vel_B[di]);
                            }
                            V_B_ptr.push_back(&p_B);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                V_B_ptr.push_back(&Y_B[si]);
                            }
                            
                            std::vector<double*> V_T_ptr;
                            V_T_ptr.push_back(&rho_T);
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                V_T_ptr.push_back(&vel_T[di]);
                            }
                            V_T_ptr.push_back(&p_T);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                V_T_ptr.push_back(&Y_T[si]);
                            }
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                V_B_ptr,
                                W_B_ptr,
                                R_y_intercell);
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                V_T_ptr,
                                W_T_ptr,
                                R_y_intercell);
                            
                            /*
                             * Convert the primitive variables into conservative variables.
                             */
                            
                            std::vector<double> m_B;
                            std::vector<double> m_T;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_B.push_back(rho_B*vel_B[di]);
                                m_T.push_back(rho_T*vel_T[di]);
                            }
                            
                            std::vector<const double*> vel_B_ptr;
                            std::vector<const double*> vel_T_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                vel_B_ptr.push_back(&vel_B[di]);
                                vel_T_ptr.push_back(&vel_T[di]);
                            }
                            
                            std::vector<const double*> Y_B_ptr;
                            std::vector<const double*> Y_T_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Y_B_ptr.push_back(&Y_B[si]);
                                Y_T_ptr.push_back(&Y_T[si]);
                            }
                            
                            double E_B = d_equation_of_state->
                                getTotalEnergyWithMassFraction(
                                    &rho_B,
                                    vel_B_ptr,
                                    &p_B,
                                    Y_B_ptr);
                            
                            double E_T = d_equation_of_state->
                                getTotalEnergyWithMassFraction(
                                    &rho_T,
                                    vel_T_ptr,
                                    &p_T,
                                    Y_T_ptr);
                            
                            /*
                             * If the WENO interpolated density or energy are negative,
                             * use constant interpolation.
                             */
                            if ((rho_B < 0) || (rho_T < 0) || (E_B < 0) || (E_T < 0))
                            {
                                // Compute the indices of bottom cell and top cell.
                                const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_T = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                rho_B = rho[idx_cell_B];
                                rho_T = rho[idx_cell_T];
                                
                                m_B[0] = rho_u[idx_cell_B];
                                m_B[1] = rho_v[idx_cell_B];
                                m_T[0] = rho_u[idx_cell_T];
                                m_T[1] = rho_v[idx_cell_T];
                                
                                E_B = E[idx_cell_B];
                                E_T = E[idx_cell_T];
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_B[si] = Y[si][idx_cell_B];
                                    Y_T[si] = Y[si][idx_cell_T];
                                }
                            }
                            
                            /*
                             * Apply the Riemann solver.
                             */
                            
                            std::vector<const double*> m_B_ptr;
                            std::vector<const double*> m_T_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_B_ptr.push_back(&m_B[di]);
                                m_T_ptr.push_back(&m_T[di]);
                            }
                            
                            std::vector<double*> F_y_midpoint_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_y_midpoint_ptr.push_back(&F_y_midpoint[ei][idx_face_y]);
                            }
                            
                            std::vector<double*> vel_y_ptr;
                            for (int vi = 0; vi < d_dim.getValue(); vi++)
                            {
                                vel_y_ptr.push_back(&(velocity_intercell->getPointer(1, vi)[idx_face_y]));
                            }
                            
                            d_Riemann_solver.computeIntercellFluxAndVelocityForFourEqnShyue(
                                F_y_midpoint_ptr,
                                vel_y_ptr,
                                &rho_B,
                                &rho_T,
                                m_B_ptr,
                                m_T_ptr,
                                &E_B,
                                &E_T,
                                Y_B_ptr,
                                Y_T_ptr,
                                Y_DIRECTION);
                        }
                    }
                    
                    // Compute the flux in the x direction.
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0] + 1; i++)
                        {
                            // Compute the indices.
                            const int idx_face_x = i +
                                j*(interior_dims[0] + 1);
                            
                            const int idx_midpoint_x = (i + 1) +
                                (j + 1)*(interior_dims[0] + 3);
                            
                            const int idx_node_L = (i - 1 + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_node_R = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            // Compute the flux.
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                convective_flux->getPointer(0, ei)[idx_face_x] = dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x + 1] +
                                    F_x_midpoint[ei][idx_midpoint_x - 1]) -
                                    3.0/10*(F_x_node[ei][idx_node_R] +
                                    F_x_node[ei][idx_node_L]) +
                                    23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
                            }
                        }
                    }
                    
                    // Compute the flux in the y direction.
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int j = 0; j < interior_dims[1] + 1; j++)
                        {
                            // Compute the indices.
                            const int idx_face_y = j +
                                i*(interior_dims[1] + 1);
                            
                            const int idx_midpoint_y = (j + 1) +
                                (i + 1)*(interior_dims[1] + 3);
                            
                            const int idx_node_B = (i + d_num_ghosts[0]) +
                                (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_node_T = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            // Compute the flux.
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                convective_flux->getPointer(1, ei)[idx_face_y] =
                                    dt*(1.0/30*(F_y_midpoint[ei][idx_midpoint_y + 1] +
                                    F_y_midpoint[ei][idx_midpoint_y - 1]) -
                                    3.0/10*(F_y_node[ei][idx_node_T] +
                                    F_y_node[ei][idx_node_B]) +
                                    23.0/15*F_y_midpoint[ei][idx_midpoint_y]);
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
                                
                                const int idx_face_x_L = (i + 1) +
                                    (j + 1)*(interior_dims[0] + 3);
                                
                                const int idx_face_x_R = (i + 2) +
                                    (j + 1)*(interior_dims[0] + 3);
                                
                                const int idx_face_y_B = (j + 1) +
                                    (i + 1)*(interior_dims[1] + 3);
                                
                                const int idx_face_y_T = (j + 2) +
                                    (i + 1)*(interior_dims[1] + 3);
                                
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
                    std::vector<double*> F_x_node;
                    std::vector<double*> F_y_node;
                    std::vector<double*> F_z_node;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
                        F_y_node.push_back(convective_flux_node[1]->getPointer(ei));
                        F_z_node.push_back(convective_flux_node[2]->getPointer(ei));
                    }
                    std::vector<double*> F_x_midpoint;
                    std::vector<double*> F_y_midpoint;
                    std::vector<double*> F_z_midpoint;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_midpoint.push_back(convective_flux_midpoint->getPointer(0, ei));
                        F_y_midpoint.push_back(convective_flux_midpoint->getPointer(1, ei));
                        F_z_midpoint.push_back(convective_flux_midpoint->getPointer(2, ei));
                    }
                    
                    // Compute the field of velocities, pressure and fluxes.
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
                                
                                std::vector<const double*> Y_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_ptr.push_back(&Y[si][idx]);
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
                                        Y_ptr);
                                
                                c[idx] = d_equation_of_state->
                                    getSoundSpeedWithMassFractionAndPressure(
                                        &rho[idx],
                                        Y_ptr,
                                        &p[idx]);
                                
                                F_x_node[0][idx] = rho_u[idx];
                                F_x_node[1][idx] = rho_u[idx]*u[idx] + p[idx];
                                F_x_node[2][idx] = rho_u[idx]*v[idx];
                                F_x_node[3][idx] = rho_u[idx]*w[idx];
                                F_x_node[4][idx] = u[idx]*(E[idx] + p[idx]);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_x_node[5 + si][idx] = u[idx]*Y[si][idx];
                                }
                                
                                F_y_node[0][idx] = rho_v[idx];
                                F_y_node[1][idx] = rho_v[idx]*u[idx];
                                F_y_node[2][idx] = rho_v[idx]*v[idx] + p[idx];
                                F_y_node[3][idx] = rho_v[idx]*w[idx];
                                F_y_node[4][idx] = v[idx]*(E[idx] + p[idx]);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_y_node[5 + si][idx] = v[idx]*Y[si][idx];
                                }
                                
                                F_z_node[0][idx] = rho_w[idx];
                                F_z_node[1][idx] = rho_w[idx]*u[idx];
                                F_z_node[2][idx] = rho_w[idx]*v[idx];
                                F_z_node[3][idx] = rho_w[idx]*w[idx] + p[idx];
                                F_z_node[4][idx] = w[idx]*(E[idx] + p[idx]);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_z_node[5 + si][idx] = w[idx]*Y[si][idx];
                                }
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
                            for (int i = -1; i < interior_dims[0] + 2; i++)
                            {
                                // Compute the indices of left cell, right cell and face of
                                // projection matrix.
                                const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_cell_R = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_face_x = (i + 1) +
                                    (j + 1)*(interior_dims[0] + 3) +
                                    (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                
                                // Get left and right quantities.
                                const double& rho_L = rho[idx_cell_L];
                                const double& rho_R = rho[idx_cell_R];
                                
                                const double& c_L = c[idx_cell_L];
                                const double& c_R = c[idx_cell_R];
                                
                                // Compute simply-averaged quantities.
                                const double rho_average = 0.5*(rho_L + rho_R);
                                const double c_average = 0.5*(c_L + c_R);
                                
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
                                            
                                        *R_x_intercell[ei][ej] = 0.0;
                                        *R_x_inv_intercell[ei][ej] = 0.0;
                                    }
                                }
                                
                                /*
                                 * Compute R_x_intercell.
                                 */
                                
                                *R_x_intercell[0][0] = 1.0/(c_average*c_average);
                                *R_x_intercell[0][1] = 1.0;
                                *R_x_intercell[0][d_num_eqn - 1] = 1.0/(c_average*c_average);
                                
                                *R_x_intercell[1][0] = -1.0/(rho_average*c_average);
                                *R_x_intercell[1][d_num_eqn - 1] = 1.0/(rho_average*c_average);
                                
                                *R_x_intercell[2][2] = 1.0;
                                
                                *R_x_intercell[3][3] = 1.0;
                                
                                *R_x_intercell[4][0] = 1.0;
                                *R_x_intercell[4][d_num_eqn - 1] = 1.0;
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    *R_x_intercell[5 + si][4 + si] = 1.0;
                                }
                                
                                /*
                                 * Compute R_x_inv_intercell.
                                 */
                                
                                *R_x_inv_intercell[0][1] = -0.5*rho_average*c_average;
                                *R_x_inv_intercell[0][4] = 0.5;
                                
                                *R_x_inv_intercell[1][0] = 1.0;
                                *R_x_inv_intercell[1][4] = -1.0/(c_average*c_average);
                                
                                *R_x_inv_intercell[2][2] = 1.0;
                                
                                *R_x_inv_intercell[3][3] = 1.0;
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    *R_x_inv_intercell[4 + si][5 + si] = 1.0;
                                }
                                
                                *R_x_inv_intercell[d_num_species + 3][1] = 0.5*rho_average*c_average;
                                *R_x_inv_intercell[d_num_species + 3][4] = 0.5;
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
                            for (int j = -1; j < interior_dims[1] + 2; j++)
                            {
                                // Compute the indices of bottom cell, top cell and face of
                                // projection matrix.
                                const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_cell_T = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_face_y = (j + 1) +
                                    (k + 1)*(interior_dims[1] + 3) +
                                    (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                
                                // Get bottom and top quantities.
                                
                                const double& rho_B = rho[idx_cell_B];
                                const double& rho_T = rho[idx_cell_T];
                                
                                const double& c_B = c[idx_cell_B];
                                const double& c_T = c[idx_cell_T];
                                
                                // Compute simply-averaged quantities.
                                const double rho_average = 0.5*(rho_B + rho_T);
                                const double c_average = 0.5*(c_B + c_T);
                                
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
                                        
                                        *R_y_intercell[ei][ej] = 0.0;
                                        *R_y_inv_intercell[ei][ej] = 0.0;
                                    }
                                }
                                
                                /*
                                 * Compute R_y_intercell.
                                 */
                                
                                *R_y_intercell[0][0] = 1.0/(c_average*c_average);
                                *R_y_intercell[0][1] = 1.0;
                                *R_y_intercell[0][d_num_eqn - 1] = 1.0/(c_average*c_average);
                                
                                *R_y_intercell[1][2] = 1.0;
                                
                                *R_y_intercell[2][0] = -1.0/(rho_average*c_average);
                                *R_y_intercell[2][d_num_eqn - 1] = 1.0/(rho_average*c_average);
                                
                                *R_y_intercell[3][3] = 1.0;
                                
                                *R_y_intercell[4][0] = 1.0;
                                *R_y_intercell[4][d_num_eqn - 1] = 1.0;
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    *R_y_intercell[5 + si][4 + si] = 1.0;
                                }
                                
                                /*
                                 * Compute R_y_inv_intercell.
                                 */
                                
                                *R_y_inv_intercell[0][2] = -0.5*rho_average*c_average;
                                *R_y_inv_intercell[0][4] = 0.5;
                                
                                *R_y_inv_intercell[1][0] = 1.0;
                                *R_y_inv_intercell[1][4] = -1.0/(c_average*c_average);
                                
                                *R_y_inv_intercell[2][1] = 1.0;
                                
                                *R_y_inv_intercell[3][3] = 1.0;
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    *R_y_inv_intercell[4 + si][5 + si] = 1.0;
                                }
                                
                                *R_y_inv_intercell[d_num_species + 3][2] = 0.5*rho_average*c_average;
                                *R_y_inv_intercell[d_num_species + 3][4] = 0.5;
                            }
                        }
                    }
                    
                    /*
                     * Compute the projection matrix and its inverse at the face normal to the
                     * z direction.
                     */
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int k = -1; k < interior_dims[2] + 2; k++)
                            {
                                // Compute the indices of back cell, front cell and face of
                                // projection matrix.
                                const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_cell_F = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_face_z = (k + 1) +
                                    (i + 1)*(interior_dims[2] + 3) +
                                    (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                
                                // Get back and front quantities.
                                
                                const double& rho_B = rho[idx_cell_B];
                                const double& rho_F = rho[idx_cell_F];
                                
                                const double& c_B = c[idx_cell_B];
                                const double& c_F = c[idx_cell_F];
                                
                                // Compute simply-averaged quantities.
                                const double rho_average = 0.5*(rho_B + rho_F);
                                const double c_average = 0.5*(c_B + c_F);
                                
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
                                        
                                        *R_z_intercell[ei][ej] = 0.0;
                                        *R_z_inv_intercell[ei][ej] = 0.0;
                                    }
                                }
                                
                                /*
                                 * Compute R_z_intercell.
                                 */
                                
                                *R_z_intercell[0][0] = 1.0/(c_average*c_average);
                                *R_z_intercell[0][1] = 1.0;
                                *R_z_intercell[0][d_num_eqn - 1] = 1.0/(c_average*c_average);
                                
                                *R_z_intercell[1][2] = 1.0;
                                
                                *R_z_intercell[2][3] = 1.0;
                                
                                *R_z_intercell[3][0] = -1.0/(rho_average*c_average);
                                *R_z_intercell[3][d_num_eqn - 1] = 1.0/(rho_average*c_average);
                                
                                *R_z_intercell[4][0] = 1.0;
                                *R_z_intercell[4][d_num_eqn - 1] = 1.0;
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    *R_z_intercell[5 + si][4 + si] = 1.0;
                                }
                                
                                /*
                                 * Compute R_z_inv_intercell.
                                 */
                                
                                *R_z_inv_intercell[0][3] = -0.5*rho_average*c_average;
                                *R_z_inv_intercell[0][4] = 0.5;
                                
                                *R_z_inv_intercell[1][0] = 1.0;
                                *R_z_inv_intercell[1][4] = -1.0/(c_average*c_average);
                                
                                *R_z_inv_intercell[2][1] = 1.0;
                                
                                *R_z_inv_intercell[3][2] = 1.0;
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    *R_z_inv_intercell[4 + si][5 + si] = 1.0;
                                }
                                
                                *R_z_inv_intercell[d_num_species + 3][3] = 0.5*rho_average*c_average;
                                *R_z_inv_intercell[d_num_species + 3][4] = 0.5;
                            }
                        }
                    }
                    
                    // Compute the mid-point flux in the x direction.
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = -1; i < interior_dims[0] + 2; i++)
                            {
                                // Compute the index of face of mid-point flux and
                                // projection matrix.
                                const int idx_face_x = (i + 1) +
                                    (j + 1)*(interior_dims[0] + 3) +
                                    (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                
                                boost::multi_array<double, 2> W_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                /*
                                 * Project primitive variables onto characteristic fields.
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
                                    
                                    std::vector<const double*> V_ptr;
                                    V_ptr.push_back(&rho[idx_cell]);
                                    V_ptr.push_back(&u[idx_cell]);
                                    V_ptr.push_back(&v[idx_cell]);
                                    V_ptr.push_back(&w[idx_cell]);
                                    V_ptr.push_back(&p[idx_cell]);
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        V_ptr.push_back(&Y[si][idx_cell]);
                                    }
                                    
                                    std::vector<double*> W_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        W_ptr.push_back(&W_array[m][ei]);
                                    }
                                    projectPrimitiveVariablesToCharacteristicFields(
                                        W_ptr,
                                        V_ptr,
                                        R_x_inv_intercell);
                                }
                                
                                /*
                                 * Do WENO interplation on characteristic variables to get W_L
                                 * and W_R
                                 */
                                
                                std::vector<double> W_L;
                                std::vector<double> W_R;
                                
                                performWENOInterpolation(W_L, W_R, W_array);
                                
                                /*
                                 * Project characteristic variables back to physical fields.
                                 */
                                
                                double rho_L;
                                double rho_R;
                                
                                std::vector<double> vel_L;
                                std::vector<double> vel_R;
                                vel_L.resize(d_dim.getValue());
                                vel_R.resize(d_dim.getValue());
                                
                                double p_L;
                                double p_R;
                                
                                std::vector<double> Y_L;
                                std::vector<double> Y_R;
                                Y_L.resize(d_num_species - 1);
                                Y_R.resize(d_num_species - 1);
                                
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
                                
                                std::vector<const double*> W_L_ptr;
                                std::vector<const double*> W_R_ptr;
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_L_ptr.push_back(&W_L[ei]);
                                    W_R_ptr.push_back(&W_R[ei]);
                                }
                                
                                std::vector<double*> V_L_ptr;
                                V_L_ptr.push_back(&rho_L);
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_L_ptr.push_back(&vel_L[di]);
                                }
                                V_L_ptr.push_back(&p_L);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_L_ptr.push_back(&Y_L[si]);
                                }
                                
                                std::vector<double*> V_R_ptr;
                                V_R_ptr.push_back(&rho_R);
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_R_ptr.push_back(&vel_R[di]);
                                }
                                V_R_ptr.push_back(&p_R);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_R_ptr.push_back(&Y_R[si]);
                                }
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_L_ptr,
                                    W_L_ptr,
                                    R_x_intercell);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_R_ptr,
                                    W_R_ptr,
                                    R_x_intercell);
                                
                                /*
                                 * Convert the primitive variables into conservative variables.
                                 */
                                
                                std::vector<double> m_L;
                                std::vector<double> m_R;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_L.push_back(rho_L*vel_L[di]);
                                    m_R.push_back(rho_R*vel_R[di]);
                                }
                                
                                std::vector<const double*> vel_L_ptr;
                                std::vector<const double*> vel_R_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    vel_L_ptr.push_back(&vel_L[di]);
                                    vel_R_ptr.push_back(&vel_R[di]);
                                }
                                
                                std::vector<const double*> Y_L_ptr;
                                std::vector<const double*> Y_R_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_L_ptr.push_back(&Y_L[si]);
                                    Y_R_ptr.push_back(&Y_R[si]);
                                }
                                
                                double E_L = d_equation_of_state->
                                    getTotalEnergyWithMassFraction(
                                        &rho_L,
                                        vel_L_ptr,
                                        &p_L,
                                        Y_L_ptr);
                                
                                double E_R = d_equation_of_state->
                                    getTotalEnergyWithMassFraction(
                                        &rho_R,
                                        vel_R_ptr,
                                        &p_R,
                                        Y_R_ptr);
                                
                                /*
                                 * If the WENO interpolated density or energy are negative,
                                 * use constant interpolation.
                                 */
                                if ((rho_L < 0) || (rho_R < 0) || (E_L < 0) || (E_R < 0))
                                {
                                    // Compute the indices of left cell and right cell.
                                    const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_R = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    rho_L = rho[idx_cell_L];
                                    rho_R = rho[idx_cell_R];
                                    
                                    m_L[0] = rho_u[idx_cell_L];
                                    m_L[1] = rho_v[idx_cell_L];
                                    m_L[2] = rho_w[idx_cell_L];
                                    m_R[0] = rho_u[idx_cell_R];
                                    m_R[1] = rho_v[idx_cell_R];
                                    m_R[2] = rho_w[idx_cell_R];
                                    
                                    E_L = E[idx_cell_L];
                                    E_R = E[idx_cell_R];
                                    
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        Y_L[si] = Y[si][idx_cell_L];
                                        Y_R[si] = Y[si][idx_cell_R];
                                    }
                                }
                                
                                /*
                                 * Apply the Riemann solver.
                                 */
                                
                                std::vector<const double*> m_L_ptr;
                                std::vector<const double*> m_R_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_L_ptr.push_back(&m_L[di]);
                                    m_R_ptr.push_back(&m_R[di]);
                                }
                                
                                std::vector<double*> F_x_midpoint_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_x_midpoint_ptr.push_back(&F_x_midpoint[ei][idx_face_x]);
                                }
                                
                                std::vector<double*> vel_x_ptr;
                                for (int vi = 0; vi < d_dim.getValue(); vi++)
                                {
                                    vel_x_ptr.push_back(&(velocity_intercell->getPointer(0, vi)[idx_face_x]));
                                }
                                
                                d_Riemann_solver.computeIntercellFluxAndVelocityForFourEqnShyue(
                                    F_x_midpoint_ptr,
                                    vel_x_ptr,
                                    &rho_L,
                                    &rho_R,
                                    m_L_ptr,
                                    m_R_ptr,
                                    &E_L,
                                    &E_R,
                                    Y_L_ptr,
                                    Y_R_ptr,
                                    X_DIRECTION);
                            }
                        }
                    }
                    
                    // Compute the mid-point flux in the y direction.
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = -1; j < interior_dims[1] + 2; j++)
                            {
                                // Compute the index of face of mid-point flux and
                                // projection matrix.
                                const int idx_face_y = (j + 1) +
                                    (k + 1)*(interior_dims[1] + 3) +
                                    (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                
                                boost::multi_array<double, 2> W_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                /*
                                 * Project primitive variables onto characteristic fields.
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
                                    
                                    std::vector<const double*> V_ptr;
                                    V_ptr.push_back(&rho[idx_cell]);
                                    V_ptr.push_back(&u[idx_cell]);
                                    V_ptr.push_back(&v[idx_cell]);
                                    V_ptr.push_back(&w[idx_cell]);
                                    V_ptr.push_back(&p[idx_cell]);
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        V_ptr.push_back(&Y[si][idx_cell]);
                                    }
                                    
                                    std::vector<double*> W_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        W_ptr.push_back(&W_array[m][ei]);
                                    }
                                    
                                    projectPrimitiveVariablesToCharacteristicFields(
                                        W_ptr,
                                        V_ptr,
                                        R_y_inv_intercell);
                                }
                                
                                /*
                                 * Do WENO interplation on characteristic variables to get W_B
                                 * and W_T
                                 */
                                
                                std::vector<double> W_B;
                                std::vector<double> W_T;
                                
                                performWENOInterpolation(W_B, W_T, W_array);
                                
                                /*
                                 * Project characteristic variables back to physical fields.
                                 */
                                
                                double rho_B;
                                double rho_T;
                                
                                std::vector<double> vel_B;
                                std::vector<double> vel_T;
                                vel_B.resize(d_dim.getValue());
                                vel_T.resize(d_dim.getValue());
                                
                                double p_B;
                                double p_T;
                                
                                std::vector<double> Y_B;
                                std::vector<double> Y_T;
                                Y_B.resize(d_num_species - 1);
                                Y_T.resize(d_num_species - 1);
                                
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
                                
                                std::vector<const double*> W_B_ptr;
                                std::vector<const double*> W_T_ptr;
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_B_ptr.push_back(&W_B[ei]);
                                    W_T_ptr.push_back(&W_T[ei]);
                                }
                                
                                std::vector<double*> V_B_ptr;
                                V_B_ptr.push_back(&rho_B);
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_B_ptr.push_back(&vel_B[di]);
                                }
                                V_B_ptr.push_back(&p_B);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_B_ptr.push_back(&Y_B[si]);
                                }
                                
                                std::vector<double*> V_T_ptr;
                                V_T_ptr.push_back(&rho_T);
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_T_ptr.push_back(&vel_T[di]);
                                }
                                V_T_ptr.push_back(&p_T);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_T_ptr.push_back(&Y_T[si]);
                                }
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_B_ptr,
                                    W_B_ptr,
                                    R_y_intercell);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_T_ptr,
                                    W_T_ptr,
                                    R_y_intercell);
                                
                                /*
                                 * Convert the primitive variables into conservative variables.
                                 */
                                
                                std::vector<double> m_B;
                                std::vector<double> m_T;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_B.push_back(rho_B*vel_B[di]);
                                    m_T.push_back(rho_T*vel_T[di]);
                                }
                                
                                std::vector<const double*> vel_B_ptr;
                                std::vector<const double*> vel_T_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    vel_B_ptr.push_back(&vel_B[di]);
                                    vel_T_ptr.push_back(&vel_T[di]);
                                }
                                
                                std::vector<const double*> Y_B_ptr;
                                std::vector<const double*> Y_T_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_B_ptr.push_back(&Y_B[si]);
                                    Y_T_ptr.push_back(&Y_T[si]);
                                }
                                
                                double E_B = d_equation_of_state->
                                    getTotalEnergyWithMassFraction(
                                        &rho_B,
                                        vel_B_ptr,
                                        &p_B,
                                        Y_B_ptr);
                                
                                double E_T = d_equation_of_state->
                                    getTotalEnergyWithMassFraction(
                                        &rho_T,
                                        vel_T_ptr,
                                        &p_T,
                                        Y_T_ptr);
                                
                                /*
                                 * If the WENO interpolated density or energy are negative,
                                 * use constant interpolation.
                                 */
                                if ((rho_B < 0) || (rho_T < 0) || (E_B < 0) || (E_T < 0))
                                {
                                    // Compute the indices of bottom cell and top cell.
                                    const int idx_cell_B = (i + d_num_ghosts[0]) +
                                        (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_T = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    rho_B = rho[idx_cell_B];
                                    rho_T = rho[idx_cell_T];
                                    
                                    m_B[0] = rho_u[idx_cell_B];
                                    m_B[1] = rho_v[idx_cell_B];
                                    m_B[2] = rho_w[idx_cell_B];
                                    m_T[0] = rho_u[idx_cell_T];
                                    m_T[1] = rho_v[idx_cell_T];
                                    m_T[2] = rho_w[idx_cell_T];
                                    
                                    E_B = E[idx_cell_B];
                                    E_T = E[idx_cell_T];
                                    
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        Y_B[si] = Y[si][idx_cell_B];
                                        Y_T[si] = Y[si][idx_cell_T];
                                    }
                                }
                                
                                /*
                                 * Apply the Riemann solver.
                                 */
                                
                                std::vector<const double*> m_B_ptr;
                                std::vector<const double*> m_T_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_B_ptr.push_back(&m_B[di]);
                                    m_T_ptr.push_back(&m_T[di]);
                                }
                                
                                std::vector<double*> F_y_midpoint_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_y_midpoint_ptr.push_back(&F_y_midpoint[ei][idx_face_y]);
                                }
                                
                                std::vector<double*> vel_y_ptr;
                                for (int vi = 0; vi < d_dim.getValue(); vi++)
                                {
                                    vel_y_ptr.push_back(&(velocity_intercell->getPointer(1, vi)[idx_face_y]));
                                }
                                
                                d_Riemann_solver.computeIntercellFluxAndVelocityForFourEqnShyue(
                                    F_y_midpoint_ptr,
                                    vel_y_ptr,
                                    &rho_B,
                                    &rho_T,
                                    m_B_ptr,
                                    m_T_ptr,
                                    &E_B,
                                    &E_T,
                                    Y_B_ptr,
                                    Y_T_ptr,
                                    Y_DIRECTION);
                            }
                        }
                    }
                    
                    // Compute the mid-point flux in the z direction.
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int k = -1; k < interior_dims[2] + 2; k++)
                            {
                                // Compute the index of face of mid-point flux and
                                // projection matrix.
                                const int idx_face_z = (k + 1) +
                                    (i + 1)*(interior_dims[2] + 3) +
                                    (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                
                                boost::multi_array<double, 2> W_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                /*
                                 * Project primitive variables onto characteristic fields.
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
                                    
                                    std::vector<const double*> V_ptr;
                                    V_ptr.push_back(&rho[idx_cell]);
                                    V_ptr.push_back(&u[idx_cell]);
                                    V_ptr.push_back(&v[idx_cell]);
                                    V_ptr.push_back(&w[idx_cell]);
                                    V_ptr.push_back(&p[idx_cell]);
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        V_ptr.push_back(&Y[si][idx_cell]);
                                    }
                                    
                                    std::vector<double*> W_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        W_ptr.push_back(&W_array[m][ei]);
                                    }
                                    
                                    projectPrimitiveVariablesToCharacteristicFields(
                                        W_ptr,
                                        V_ptr,
                                        R_z_inv_intercell);
                                }
                                
                                /*
                                 * Do WENO interplation on characteristic variables to get W_B
                                 * and W_F
                                 */
                                
                                std::vector<double> W_B;
                                std::vector<double> W_F;
                                
                                performWENOInterpolation(W_B, W_F, W_array);
                                
                                /*
                                 * Project characteristic variables back to physical fields.
                                 */
                                
                                double rho_B;
                                double rho_F;
                                
                                std::vector<double> vel_B;
                                std::vector<double> vel_F;
                                vel_B.resize(d_dim.getValue());
                                vel_F.resize(d_dim.getValue());
                                
                                double p_B;
                                double p_F;
                                
                                std::vector<double> Y_B;
                                std::vector<double> Y_F;
                                Y_B.resize(d_num_species - 1);
                                Y_F.resize(d_num_species - 1);
                                
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
                                
                                std::vector<const double*> W_B_ptr;
                                std::vector<const double*> W_F_ptr;
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_B_ptr.push_back(&W_B[ei]);
                                    W_F_ptr.push_back(&W_F[ei]);
                                }
                                
                                std::vector<double*> V_B_ptr;
                                V_B_ptr.push_back(&rho_B);
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_B_ptr.push_back(&vel_B[di]);
                                }
                                V_B_ptr.push_back(&p_B);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_B_ptr.push_back(&Y_B[si]);
                                }
                                
                                std::vector<double*> V_F_ptr;
                                V_F_ptr.push_back(&rho_F);
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_F_ptr.push_back(&vel_F[di]);
                                }
                                V_F_ptr.push_back(&p_F);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_F_ptr.push_back(&Y_F[si]);
                                }
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_B_ptr,
                                    W_B_ptr,
                                    R_z_intercell);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_F_ptr,
                                    W_F_ptr,
                                    R_z_intercell);
                                
                                /*
                                 * Convert the primitive variables into conservative variables.
                                 */
                                
                                
                                std::vector<double> m_B;
                                std::vector<double> m_F;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_B.push_back(rho_B*vel_B[di]);
                                    m_F.push_back(rho_F*vel_F[di]);
                                }
                                
                                std::vector<const double*> vel_B_ptr;
                                std::vector<const double*> vel_F_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    vel_B_ptr.push_back(&vel_B[di]);
                                    vel_F_ptr.push_back(&vel_F[di]);
                                }
                                
                                std::vector<const double*> Y_B_ptr;
                                std::vector<const double*> Y_F_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Y_B_ptr.push_back(&Y_B[si]);
                                    Y_F_ptr.push_back(&Y_F[si]);
                                }
                                
                                double E_B = d_equation_of_state->
                                    getTotalEnergyWithMassFraction(
                                        &rho_B,
                                        vel_B_ptr,
                                        &p_B,
                                        Y_B_ptr);
                                
                                double E_F = d_equation_of_state->
                                    getTotalEnergyWithMassFraction(
                                        &rho_F,
                                        vel_F_ptr,
                                        &p_F,
                                        Y_F_ptr);
                                
                                /*
                                 * If the WENO interpolated density or energy are negative,
                                 * use constant interpolation.
                                 */
                                if ((rho_B < 0) || (rho_F < 0) || (E_B < 0) || (E_F < 0))
                                {
                                    // Compute the indices of back cell and front cell.
                                    const int idx_cell_B = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_F = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                        
                                    rho_B = rho[idx_cell_B];
                                    rho_F = rho[idx_cell_F];
                                    
                                    m_B[0] = rho_u[idx_cell_B];
                                    m_B[1] = rho_v[idx_cell_B];
                                    m_B[2] = rho_w[idx_cell_B];
                                    m_F[0] = rho_u[idx_cell_F];
                                    m_F[1] = rho_v[idx_cell_F];
                                    m_F[2] = rho_w[idx_cell_F];
                                    
                                    E_B = E[idx_cell_B];
                                    E_F = E[idx_cell_F];
                                    
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        Y_B[si] = Y[si][idx_cell_B];
                                        Y_F[si] = Y[si][idx_cell_F];
                                    }
                                }
                                
                                /*
                                 * Apply the Riemann solver.
                                 */
                                
                                std::vector<const double*> m_B_ptr;
                                std::vector<const double*> m_F_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_B_ptr.push_back(&m_B[di]);
                                    m_F_ptr.push_back(&m_F[di]);
                                }
                                
                                std::vector<double*> F_z_midpoint_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_z_midpoint_ptr.push_back(&F_z_midpoint[ei][idx_face_z]);
                                }
                                
                                std::vector<double*> vel_z_ptr;
                                for (int vi = 0; vi < d_dim.getValue(); vi++)
                                {
                                    vel_z_ptr.push_back(&(velocity_intercell->getPointer(2, vi)[idx_face_z]));
                                }
                                
                                d_Riemann_solver.computeIntercellFluxAndVelocityForFourEqnShyue(
                                    F_z_midpoint_ptr,
                                    vel_z_ptr,
                                    &rho_B,
                                    &rho_F,
                                    m_B_ptr,
                                    m_F_ptr,
                                    &E_B,
                                    &E_F,
                                    Y_B_ptr,
                                    Y_F_ptr,
                                    Z_DIRECTION);
                            }
                        }
                    }
                    
                    // Compute the flux in the x direction.
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0] + 1; i++)
                            {
                                // Compute the indices.
                                const int idx_face_x = i +
                                    j*(interior_dims[0] + 1) +
                                    k*(interior_dims[0] + 1)*interior_dims[1];
                                
                                const int idx_midpoint_x = (i + 1) +
                                    (j + 1)*(interior_dims[0] + 3) +
                                    (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                
                                const int idx_node_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_node_R = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                // Compute the flux.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    convective_flux->getPointer(0, ei)[idx_face_x] = dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x + 1] +
                                        F_x_midpoint[ei][idx_midpoint_x - 1]) -
                                        3.0/10*(F_x_node[ei][idx_node_R] +
                                        F_x_node[ei][idx_node_L]) +
                                        23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
                                }
                            }
                        }
                    }
                    
                    // Compute the flux in the y direction.
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = 0; j < interior_dims[1] + 1; j++)
                            {
                                // Compute the indices.
                                const int idx_face_y = j +
                                    k*(interior_dims[1] + 1) +
                                    i*(interior_dims[1] + 1)*interior_dims[2];
                                
                                const int idx_midpoint_y = (j + 1) +
                                    (k + 1)*(interior_dims[1] + 3) +
                                    (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                
                                const int idx_node_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_node_T = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                // Compute the flux.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    convective_flux->getPointer(1, ei)[idx_face_y] =
                                        dt*(1.0/30*(F_y_midpoint[ei][idx_midpoint_y + 1] +
                                        F_y_midpoint[ei][idx_midpoint_y - 1]) -
                                        3.0/10*(F_y_node[ei][idx_node_T] +
                                        F_y_node[ei][idx_node_B]) +
                                        23.0/15*F_y_midpoint[ei][idx_midpoint_y]);
                                }
                            }
                        }
                    }
                    
                    // Compute the flux in the z direction.
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int k = 0; k < interior_dims[2] + 1; k++)
                            {
                                // Compute the indices.
                                const int idx_face_z = k +
                                    i*(interior_dims[2] + 1) +
                                    j*(interior_dims[2] + 1)*interior_dims[0];
                                
                                const int idx_midpoint_z = (k + 1) +
                                    (i + 1)*(interior_dims[2] + 3) +
                                    (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                
                                const int idx_node_B = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_node_F = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                // Compute the flux.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    convective_flux->getPointer(2, ei)[idx_face_z] =
                                        dt*(1.0/30*(F_z_midpoint[ei][idx_midpoint_z + 1] +
                                        F_z_midpoint[ei][idx_midpoint_z - 1]) -
                                        3.0/10*(F_z_node[ei][idx_node_F] +
                                        F_z_node[ei][idx_node_B]) +
                                        23.0/15*F_z_midpoint[ei][idx_midpoint_z]);
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
                                    
                                    const int idx_face_x_L = (i + 1) +
                                        (j + 1)*(interior_dims[0] + 3) +
                                        (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                    
                                    const int idx_face_x_R = (i + 2) +
                                        (j + 1)*(interior_dims[0] + 3) +
                                        (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                    
                                    const int idx_face_y_B = (j + 1) +
                                        (k + 1)*(interior_dims[1] + 3) +
                                        (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                    
                                    const int idx_face_y_T = (j + 2) +
                                        (k + 1)*(interior_dims[1] + 3) +
                                        (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                    
                                    const int idx_face_z_B = (k + 1) +
                                        (i + 1)*(interior_dims[2] + 3) +
                                        (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                    
                                    const int idx_face_z_F = (k + 2) +
                                        (i + 1)*(interior_dims[2] + 3) +
                                        (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                    
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
                    new pdat::FaceData<double>(interior_box, d_dim.getValue(), hier::IntVector::getOne(d_dim)));
                
                boost::shared_ptr<pdat::CellData<double> > density(
                    new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                
                boost::shared_ptr<pdat::CellData<double> > velocity(
                    new pdat::CellData<double>(interior_box, d_dim.getValue(), d_num_ghosts));
                
                boost::shared_ptr<pdat::CellData<double> > pressure(
                    new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                
                boost::shared_ptr<pdat::CellData<double> > sound_speed(
                    new pdat::CellData<double>(interior_box, 1, d_num_ghosts));
                
                std::vector<boost::shared_ptr<pdat::CellData<double> > > convective_flux_node;
                for (int di = 0; di < d_dim.getValue(); di++)
                {
                    convective_flux_node.push_back(boost::make_shared<pdat::CellData<double> >(
                        interior_box, d_num_eqn, d_num_ghosts));
                }
                
                boost::shared_ptr<pdat::FaceData<double> > convective_flux_midpoint(
                    new pdat::FaceData<double>(interior_box, d_num_eqn, hier::IntVector::getOne(d_dim)));
                
                boost::shared_ptr<pdat::FaceData<double> > projection_matrix(
                    new pdat::FaceData<double>(
                        interior_box,
                        d_num_eqn*d_num_eqn,
                        hier::IntVector::getOne(d_dim)));
                
                boost::shared_ptr<pdat::FaceData<double> > projection_matrix_inv(
                    new pdat::FaceData<double>(
                        interior_box,
                        d_num_eqn*d_num_eqn,
                        hier::IntVector::getOne(d_dim)));
                
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
                    std::vector<double*> F_x_node;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
                    }
                    std::vector<double*> F_x_midpoint;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_midpoint.push_back(convective_flux_midpoint->getPointer(0, ei));
                    }
                    
                    // Compute the field of total density, velocities, pressure and fluxes.
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
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            F_x_node[si][idx] = u[idx]*Z_rho[si][idx];
                        }
                        F_x_node[d_num_species][idx] = rho_u[idx]*u[idx] + p[idx];
                        F_x_node[d_num_species + 1][idx] = u[idx]*(E[idx] + p[idx]);
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            F_x_node[d_num_species + 2 + si][idx] = u[idx]*Z[si][idx];
                        }
                    }
                    
                    /*
                     * Compute the projection matrix and its inverse at the face normal to the
                     * x direction.
                     */
                    for (int i = -1; i < interior_dims[0] + 2; i++)
                    {
                        // Compute the indices of left cell, right cell and face of
                        // projection matrix.
                        const int idx_cell_L = i - 1 + d_num_ghosts[0];
                        
                        const int idx_cell_R = i + d_num_ghosts[0];
                        
                        const int idx_face_x = i + 1;
                        
                        // Get left and right quantities.
                        const double& rho_L = rho[idx_cell_L];
                        const double& rho_R = rho[idx_cell_R];
                        
                        const double& c_L = c[idx_cell_L];
                        const double& c_R = c[idx_cell_R];
                        
                        // Compute simply-averaged quantities.
                        const double rho_average = 0.5*(rho_L + rho_R);
                        const double c_average = 0.5*(c_L + c_R);
                        
                        std::vector<double> Z_rho_average;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho_average.push_back(0.5*(Z_rho[si][idx_cell_L] +
                                Z_rho[si][idx_cell_R]));
                        }
                        
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
                                    
                                *R_x_intercell[ei][ej] = 0.0;
                                *R_x_inv_intercell[ei][ej] = 0.0;
                            }
                        }
                        
                        /*
                         * Compute R_x_intercell.
                         */
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            *R_x_intercell[si][0] = -0.5*Z_rho_average[si]/c_average;
                            *R_x_intercell[si][si + 1] = 1.0;
                            *R_x_intercell[si][d_num_eqn - 1] = 0.5*Z_rho_average[si]/c_average;
                        }
                        
                        *R_x_intercell[d_num_species][0] = 0.5;
                        *R_x_intercell[d_num_species][d_num_eqn - 1] = 0.5;
                        
                        *R_x_intercell[d_num_species + 1][0] = -0.5*rho_average*c_average;
                        *R_x_intercell[d_num_species + 1][d_num_eqn - 1] = 0.5*rho_average*c_average;
                        
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            *R_x_intercell[d_num_species + 2 + si][d_num_species + 1 + si] = 1.0;
                        }
                        
                        /*
                         * Compute R_x_inv_intercell.
                         */
                        
                        *R_x_inv_intercell[0][d_num_species] = 1.0;
                        *R_x_inv_intercell[0][d_num_species + 1] = -1.0/(rho_average*c_average);
                        
                        for (int si = 0; si < d_num_species; si++)
                        {
                            *R_x_inv_intercell[1 + si][si] = 1.0;
                            *R_x_inv_intercell[1 + si][d_num_species + 1] = -Z_rho_average[si]/
                                (rho_average*c_average*c_average);
                        }
                        
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            *R_x_inv_intercell[d_num_species + 1 + si][d_num_species + 2 + si] = 1.0;
                        }
                        
                        *R_x_inv_intercell[2*d_num_species][d_num_species] = 1.0;
                        *R_x_inv_intercell[2*d_num_species][d_num_species + 1] = 1.0/(rho_average*c_average);
                    }
                    
                    // Compute the mid-point flux in the x direction.
                    for (int i = -1; i < interior_dims[0] + 2; i++)
                    {
                        // Compute the index of face of mid-point flux and
                        // projection matrix.
                        const int idx_face_x = i + 1;
                        
                        boost::multi_array<double, 2> W_array(
                            boost::extents[6][d_num_eqn]);
                        
                        /*
                         * Project primitive variables onto characteristic fields.
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
                            
                            std::vector<const double*> V_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                V_ptr.push_back(&Z_rho[si][idx_cell]);
                            }
                            V_ptr.push_back(&u[idx_cell]);
                            V_ptr.push_back(&p[idx_cell]);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                V_ptr.push_back(&Z[si][idx_cell]);
                            }
                            
                            std::vector<double*> W_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                W_ptr.push_back(&W_array[m][ei]);
                            }
                            projectPrimitiveVariablesToCharacteristicFields(
                                W_ptr,
                                V_ptr,
                                R_x_inv_intercell);
                        }
                        
                        /*
                         * Do WENO interplation on characteristic variables to get W_L
                         * and W_R
                         */
                        
                        std::vector<double> W_L;
                        std::vector<double> W_R;
                        
                        performWENOInterpolation(W_L, W_R, W_array);
                        
                        /*
                         * Project characteristic variables back to physical fields.
                         */
                        
                        std::vector<double> Z_rho_L;
                        std::vector<double> Z_rho_R;
                        Z_rho_L.resize(d_num_species);
                        Z_rho_R.resize(d_num_species);
                        
                        std::vector<double> vel_L;
                        std::vector<double> vel_R;
                        vel_L.resize(d_dim.getValue());
                        vel_R.resize(d_dim.getValue());
                        
                        double p_L;
                        double p_R;
                        
                        std::vector<double> Z_L;
                        std::vector<double> Z_R;
                        Z_L.resize(d_num_species - 1);
                        Z_R.resize(d_num_species - 1);
                        
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
                        
                        std::vector<const double*> W_L_ptr;
                        std::vector<const double*> W_R_ptr;
                        
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            W_L_ptr.push_back(&W_L[ei]);
                            W_R_ptr.push_back(&W_R[ei]);
                        }
                        
                        std::vector<double*> V_L_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            V_L_ptr.push_back(&Z_rho_L[si]);
                        }
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            V_L_ptr.push_back(&vel_L[di]);
                        }
                        V_L_ptr.push_back(&p_L);
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            V_L_ptr.push_back(&Z_L[si]);
                        }
                        
                        std::vector<double*> V_R_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            V_R_ptr.push_back(&Z_rho_R[si]);
                        }
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            V_R_ptr.push_back(&vel_R[di]);
                        }
                        V_R_ptr.push_back(&p_R);
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            V_R_ptr.push_back(&Z_R[si]);
                        }
                        
                        projectCharacteristicVariablesToPhysicalFields(
                            V_L_ptr,
                            W_L_ptr,
                            R_x_intercell);
                        
                        projectCharacteristicVariablesToPhysicalFields(
                            V_R_ptr,
                            W_R_ptr,
                            R_x_intercell);
                        
                        /*
                         * Convert the primitive variables into conservative variables.
                         */
                        
                        std::vector<const double*> Z_rho_L_ptr;
                        std::vector<const double*> Z_rho_R_ptr;
                        for (int si = 0; si < d_num_species; si++)
                        {
                            Z_rho_L_ptr.push_back(&Z_rho_L[si]);
                            Z_rho_R_ptr.push_back(&Z_rho_R[si]);
                        }
                        
                        double rho_L = d_equation_of_state->
                            getTotalDensity(
                                Z_rho_L_ptr);
                        
                        double rho_R = d_equation_of_state->
                            getTotalDensity(
                                Z_rho_R_ptr);
                        
                        std::vector<double> m_L;
                        std::vector<double> m_R;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            m_L.push_back(rho_L*vel_L[di]);
                            m_R.push_back(rho_R*vel_R[di]);
                        }
                        
                        std::vector<const double*> vel_L_ptr;
                        std::vector<const double*> vel_R_ptr;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            vel_L_ptr.push_back(&vel_L[di]);
                            vel_R_ptr.push_back(&vel_R[di]);
                        }
                        
                        std::vector<const double*> Z_L_ptr;
                        std::vector<const double*> Z_R_ptr;
                        for (int si = 0; si < d_num_species - 1; si++)
                        {
                            Z_L_ptr.push_back(&Z_L[si]);
                            Z_R_ptr.push_back(&Z_R[si]);
                        }
                        
                        double E_L = d_equation_of_state->
                            getTotalEnergyWithVolumeFraction(
                                &rho_L,
                                vel_L_ptr,
                                &p_L,
                                Z_L_ptr);
                        
                        double E_R = d_equation_of_state->
                            getTotalEnergyWithVolumeFraction(
                                &rho_R,
                                vel_R_ptr,
                                &p_R,
                                Z_R_ptr);
                        
                        /*
                         * If the WENO interpolated density or energy are negative,
                         * use constant interpolation.
                         */
                        if ((rho_L < 0) || (rho_R < 0) || (E_L < 0) || (E_R < 0))
                        {
                            // Compute the indices of left cell and right cell.
                            const int idx_cell_L = i - 1 + d_num_ghosts[0];
                            const int idx_cell_R = i + d_num_ghosts[0];
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_L[si] = Z_rho[si][idx_cell_L];
                                Z_rho_R[si] = Z_rho[si][idx_cell_R];
                            }
                            
                            m_L[0] = rho_u[idx_cell_L];
                            m_R[0] = rho_u[idx_cell_R];
                            
                            E_L = E[idx_cell_L];
                            E_R = E[idx_cell_R];
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Z_L[si] = Z[si][idx_cell_L];
                                Z_R[si] = Z[si][idx_cell_R];
                            }
                        }
                        
                        /*
                         * Apply the Riemann solver.
                         */
                        
                        std::vector<const double*> m_L_ptr;
                        std::vector<const double*> m_R_ptr;
                        for (int di = 0; di < d_dim.getValue(); di++)
                        {
                            m_L_ptr.push_back(&m_L[di]);
                            m_R_ptr.push_back(&m_R[di]);
                        }
                        
                        std::vector<double*> F_x_midpoint_ptr;
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            F_x_midpoint_ptr.push_back(&F_x_midpoint[ei][idx_face_x]);
                        }
                        
                        std::vector<double*> vel_x_ptr;
                        for (int vi = 0; vi < d_dim.getValue(); vi++)
                        {
                            vel_x_ptr.push_back(&(velocity_intercell->getPointer(0, vi)[idx_face_x]));
                        }
                        
                        d_Riemann_solver.computeIntercellFluxAndVelocityForFiveEqnAllaire(
                            F_x_midpoint_ptr,
                            vel_x_ptr,
                            Z_rho_L_ptr,
                            Z_rho_R_ptr,
                            m_L_ptr,
                            m_R_ptr,
                            &E_L,
                            &E_R,
                            Z_L_ptr,
                            Z_R_ptr,
                            X_DIRECTION);
                    }
                    
                    // Compute the flux in the x direction.
                    for (int i = 0; i < interior_dims[0] + 1; i++)
                    {
                        // Compute the indices.
                        const int idx_face_x = i;
                        
                        const int idx_midpoint_x = i + 1;
                        
                        const int idx_node_L = i - 1 + d_num_ghosts[0];
                        
                        const int idx_node_R = i + d_num_ghosts[0];
                        
                        // Compute the flux.
                        for (int ei = 0; ei < d_num_eqn; ei++)
                        {
                            convective_flux->getPointer(0, ei)[idx_face_x] = dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x + 1] +
                                F_x_midpoint[ei][idx_midpoint_x - 1]) -
                                3.0/10*(F_x_node[ei][idx_node_R] +
                                F_x_node[ei][idx_node_L]) +
                                23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
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
                            
                            const int idx_face_x_L = i + 1;
                            
                            const int idx_face_x_R = i + 2;
                            
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
                    std::vector<double*> F_x_node;
                    std::vector<double*> F_y_node;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
                        F_y_node.push_back(convective_flux_node[1]->getPointer(ei));
                    }
                    std::vector<double*> F_x_midpoint;
                    std::vector<double*> F_y_midpoint;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_midpoint.push_back(convective_flux_midpoint->getPointer(0, ei));
                        F_y_midpoint.push_back(convective_flux_midpoint->getPointer(1, ei));
                    }
                    
                    // Compute the field of total density, velocities, pressure and fluxes.
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
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                F_x_node[si][idx] = u[idx]*Z_rho[si][idx];
                            }
                            F_x_node[d_num_species][idx] = rho_u[idx]*u[idx] + p[idx];
                            F_x_node[d_num_species + 1][idx] = rho_u[idx]*v[idx];
                            F_x_node[d_num_species + 2][idx] = u[idx]*(E[idx] + p[idx]);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_x_node[d_num_species + 3 + si][idx] = u[idx]*Z[si][idx];
                            }
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                F_y_node[si][idx] = v[idx]*Z_rho[si][idx];
                            }
                            F_y_node[d_num_species][idx] = rho_v[idx]*u[idx];
                            F_y_node[d_num_species + 1][idx] = rho_v[idx]*v[idx] + p[idx];
                            F_y_node[d_num_species + 2][idx] = v[idx]*(E[idx] + p[idx]);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                F_y_node[d_num_species + 3 + si][idx] = v[idx]*Z[si][idx];
                            }
                        }
                    }
                    
                    /*
                     * Compute the projection matrix and its inverse at the face normal to the
                     * x direction.
                     */
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = -1; i < interior_dims[0] + 2; i++)
                        {
                            // Compute the indices of left cell, right cell and face of
                            // projection matrix.
                            const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_cell_R = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_face_x = (i + 1) +
                                (j + 1)*(interior_dims[0] + 3);
                            
                            // Get left and right quantities.
                            const double& rho_L = rho[idx_cell_L];
                            const double& rho_R = rho[idx_cell_R];
                            
                            const double& c_L = c[idx_cell_L];
                            const double& c_R = c[idx_cell_R];
                            
                            // Compute simply-averaged quantities.
                            const double rho_average = 0.5*(rho_L + rho_R);
                            const double c_average = 0.5*(c_L + c_R);
                            
                            std::vector<double> Z_rho_average;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_average.push_back(0.5*(Z_rho[si][idx_cell_L] +
                                    Z_rho[si][idx_cell_R]));
                            }
                            
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
                                        
                                    *R_x_intercell[ei][ej] = 0.0;
                                    *R_x_inv_intercell[ei][ej] = 0.0;
                                }
                            }
                            
                            /*
                             * Compute R_x_intercell.
                             */
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                *R_x_intercell[si][0] = -0.5*Z_rho_average[si]/c_average;
                                *R_x_intercell[si][si + 1] = 1.0;
                                *R_x_intercell[si][d_num_eqn - 1] = 0.5*Z_rho_average[si]/c_average;
                            }
                            
                            *R_x_intercell[d_num_species][0] = 0.5;
                            *R_x_intercell[d_num_species][d_num_eqn - 1] = 0.5;
                            
                            *R_x_intercell[d_num_species + 1][d_num_species + 1] = 1.0;
                            
                            *R_x_intercell[d_num_species + 2][0] = -0.5*rho_average*c_average;
                            *R_x_intercell[d_num_species + 2][d_num_eqn - 1] = 0.5*rho_average*c_average;
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                *R_x_intercell[d_num_species + 3 + si][d_num_species + 2 + si] = 1.0;
                            }
                            
                            /*
                             * Compute R_x_inv_intercell.
                             */
                            
                            *R_x_inv_intercell[0][d_num_species] = 1.0;
                            *R_x_inv_intercell[0][d_num_species + 2] = -1.0/(rho_average*c_average);
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                *R_x_inv_intercell[1 + si][si] = 1.0;
                                *R_x_inv_intercell[1 + si][d_num_species + 2] = -Z_rho_average[si]/
                                    (rho_average*c_average*c_average);
                            }
                            
                            *R_x_inv_intercell[d_num_species + 1][d_num_species + 1] = 1.0;
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                *R_x_inv_intercell[d_num_species + 2 + si][d_num_species + 3 + si] = 1.0;
                            }
                            
                            *R_x_inv_intercell[2*d_num_species + 1][d_num_species] = 1.0;
                            *R_x_inv_intercell[2*d_num_species + 1][d_num_species + 2] = 1.0/(rho_average*c_average);
                        }
                    }
                    
                    /*
                     * Compute the projection matrix and its inverse at the face normal to the
                     * y direction.
                     */
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int j = -1; j < interior_dims[1] + 2; j++)
                        {
                            // Compute the indices of bottom cell, top cell and face of
                            // projection matrix.
                            const int idx_cell_B = (i + d_num_ghosts[0]) +
                                (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_cell_T = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_face_y = (j + 1) +
                                (i + 1)*(interior_dims[1] + 3);
                            
                            // Get bottom and top quantities.
                            
                            const double& rho_B = rho[idx_cell_B];
                            const double& rho_T = rho[idx_cell_T];
                            
                            const double& c_B = c[idx_cell_B];
                            const double& c_T = c[idx_cell_T];
                            
                            // Compute simply-averaged quantities.
                            const double rho_average = 0.5*(rho_B + rho_T);
                            const double c_average = 0.5*(c_B + c_T);
                            
                            std::vector<double> Z_rho_average;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_average.push_back(0.5*(Z_rho[si][idx_cell_B] +
                                    Z_rho[si][idx_cell_T]));
                            }
                            
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
                                    
                                    *R_y_intercell[ei][ej] = 0.0;
                                    *R_y_inv_intercell[ei][ej] = 0.0;
                                }
                            }
                            
                            /*
                             * Compute R_y_intercell.
                             */
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                *R_y_intercell[si][0] = -0.5*Z_rho_average[si]/c_average;
                                *R_y_intercell[si][si + 1] = 1.0;
                                *R_y_intercell[si][d_num_eqn - 1] = 0.5*Z_rho_average[si]/c_average;
                            }
                            
                            *R_y_intercell[d_num_species][d_num_species + 1] = 1.0;
                            
                            *R_y_intercell[d_num_species + 1][0] = 0.5;
                            *R_y_intercell[d_num_species + 1][d_num_eqn - 1] = 0.5;
                            
                            *R_y_intercell[d_num_species + 2][0] = -0.5*rho_average*c_average;
                            *R_y_intercell[d_num_species + 2][d_num_eqn - 1] = 0.5*rho_average*c_average;
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                *R_y_intercell[d_num_species + 3 + si][d_num_species + 2 + si] = 1.0;
                            }
                            
                            /*
                             * Compute R_y_inv_intercell.
                             */
                            
                            *R_y_inv_intercell[0][d_num_species + 1] = 1.0;
                            *R_y_inv_intercell[0][d_num_species + 2] = -1.0/(rho_average*c_average);
                            
                            for (int si = 0; si < d_num_species; si++)
                            {
                                *R_y_inv_intercell[1 + si][si] = 1.0;
                                *R_y_inv_intercell[1 + si][d_num_species + 2] = -Z_rho_average[si]/
                                    (rho_average*c_average*c_average);
                            }
                            
                            *R_y_inv_intercell[d_num_species + 1][d_num_species] = 1.0;
                            
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                *R_y_inv_intercell[d_num_species + 2 + si][d_num_species + 3 + si] = 1.0;
                            }
                            
                            *R_y_inv_intercell[2*d_num_species + 1][d_num_species + 1] = 1.0;
                            *R_y_inv_intercell[2*d_num_species + 1][d_num_species + 2] = 1.0/(rho_average*c_average);
                        }
                    }
                    
                    // Compute the mid-point flux in the x direction.
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = -1; i < interior_dims[0] + 2; i++)
                        {
                            // Compute the index of face of mid-point flux and
                            // projection matrix.
                            const int idx_face_x = (i + 1) +
                                (j + 1)*(interior_dims[0] + 3);
                            
                            boost::multi_array<double, 2> W_array(
                                boost::extents[6][d_num_eqn]);
                            
                            /*
                             * Project primitive variables onto characteristic fields.
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
                                
                                std::vector<const double*> V_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ptr.push_back(&Z_rho[si][idx_cell]);
                                }
                                V_ptr.push_back(&u[idx_cell]);
                                V_ptr.push_back(&v[idx_cell]);
                                V_ptr.push_back(&p[idx_cell]);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_ptr.push_back(&Z[si][idx_cell]);
                                }
                                
                                std::vector<double*> W_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_ptr.push_back(&W_array[m][ei]);
                                }
                                projectPrimitiveVariablesToCharacteristicFields(
                                    W_ptr,
                                    V_ptr,
                                    R_x_inv_intercell);
                            }
                            
                            /*
                             * Do WENO interplation on characteristic variables to get W_L
                             * and W_R
                             */
                            
                            std::vector<double> W_L;
                            std::vector<double> W_R;
                            
                            performWENOInterpolation(W_L, W_R, W_array);
                            
                            /*
                             * Project characteristic variables back to physical fields.
                             */
                            
                            std::vector<double> Z_rho_L;
                            std::vector<double> Z_rho_R;
                            Z_rho_L.resize(d_num_species);
                            Z_rho_R.resize(d_num_species);
                            
                            std::vector<double> vel_L;
                            std::vector<double> vel_R;
                            vel_L.resize(d_dim.getValue());
                            vel_R.resize(d_dim.getValue());
                            
                            double p_L;
                            double p_R;
                            
                            std::vector<double> Z_L;
                            std::vector<double> Z_R;
                            Z_L.resize(d_num_species - 1);
                            Z_R.resize(d_num_species - 1);
                            
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
                            
                            std::vector<const double*> W_L_ptr;
                            std::vector<const double*> W_R_ptr;
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                W_L_ptr.push_back(&W_L[ei]);
                                W_R_ptr.push_back(&W_R[ei]);
                            }
                            
                            std::vector<double*> V_L_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                V_L_ptr.push_back(&Z_rho_L[si]);
                            }
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                V_L_ptr.push_back(&vel_L[di]);
                            }
                            V_L_ptr.push_back(&p_L);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                V_L_ptr.push_back(&Z_L[si]);
                            }
                            
                            std::vector<double*> V_R_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                V_R_ptr.push_back(&Z_rho_R[si]);
                            }
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                V_R_ptr.push_back(&vel_R[di]);
                            }
                            V_R_ptr.push_back(&p_R);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                V_R_ptr.push_back(&Z_R[si]);
                            }
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                V_L_ptr,
                                W_L_ptr,
                                R_x_intercell);
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                V_R_ptr,
                                W_R_ptr,
                                R_x_intercell);
                            
                            /*
                             * Convert the primitive variables into conservative variables.
                             */
                            
                            std::vector<const double*> Z_rho_L_ptr;
                            std::vector<const double*> Z_rho_R_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_L_ptr.push_back(&Z_rho_L[si]);
                                Z_rho_R_ptr.push_back(&Z_rho_R[si]);
                            }
                            
                            double rho_L = d_equation_of_state->
                                getTotalDensity(
                                    Z_rho_L_ptr);
                            
                            double rho_R = d_equation_of_state->
                                getTotalDensity(
                                    Z_rho_R_ptr);
                            
                            std::vector<double> m_L;
                            std::vector<double> m_R;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_L.push_back(rho_L*vel_L[di]);
                                m_R.push_back(rho_R*vel_R[di]);
                            }
                            
                            std::vector<const double*> vel_L_ptr;
                            std::vector<const double*> vel_R_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                vel_L_ptr.push_back(&vel_L[di]);
                                vel_R_ptr.push_back(&vel_R[di]);
                            }
                            
                            std::vector<const double*> Z_L_ptr;
                            std::vector<const double*> Z_R_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Z_L_ptr.push_back(&Z_L[si]);
                                Z_R_ptr.push_back(&Z_R[si]);
                            }
                            
                            double E_L = d_equation_of_state->
                                getTotalEnergyWithVolumeFraction(
                                    &rho_L,
                                    vel_L_ptr,
                                    &p_L,
                                    Z_L_ptr);
                            
                            double E_R = d_equation_of_state->
                                getTotalEnergyWithVolumeFraction(
                                    &rho_R,
                                    vel_R_ptr,
                                    &p_R,
                                    Z_R_ptr);
                            
                            /*
                             * If the WENO interpolated density or energy are negative,
                             * use constant interpolation.
                             */
                            if ((rho_L < 0) || (rho_R < 0) || (E_L < 0) || (E_R < 0))
                            {
                                // Compute the indices of left cell and right cell.
                                const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_R = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_L[si] = Z_rho[si][idx_cell_L];
                                    Z_rho_R[si] = Z_rho[si][idx_cell_R];
                                }
                                
                                m_L[0] = rho_u[idx_cell_L];
                                m_L[1] = rho_v[idx_cell_L];
                                m_R[0] = rho_u[idx_cell_R];
                                m_R[1] = rho_v[idx_cell_R];
                                
                                E_L = E[idx_cell_L];
                                E_R = E[idx_cell_R];
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_L[si] = Z[si][idx_cell_L];
                                    Z_R[si] = Z[si][idx_cell_R];
                                }
                            }
                            
                            /*
                             * Apply the Riemann solver.
                             */
                            
                            std::vector<const double*> m_L_ptr;
                            std::vector<const double*> m_R_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_L_ptr.push_back(&m_L[di]);
                                m_R_ptr.push_back(&m_R[di]);
                            }
                            
                            std::vector<double*> F_x_midpoint_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_x_midpoint_ptr.push_back(&F_x_midpoint[ei][idx_face_x]);
                            }
                            
                            std::vector<double*> vel_x_ptr;
                            for (int vi = 0; vi < d_dim.getValue(); vi++)
                            {
                                vel_x_ptr.push_back(&(velocity_intercell->getPointer(0, vi)[idx_face_x]));
                            }
                            
                            d_Riemann_solver.computeIntercellFluxAndVelocityForFiveEqnAllaire(
                                F_x_midpoint_ptr,
                                vel_x_ptr,
                                Z_rho_L_ptr,
                                Z_rho_R_ptr,
                                m_L_ptr,
                                m_R_ptr,
                                &E_L,
                                &E_R,
                                Z_L_ptr,
                                Z_R_ptr,
                                X_DIRECTION);
                        }
                    }
                    
                    // Compute the mid-point flux in the y direction.
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int j = -1; j < interior_dims[1] + 2; j++)
                        {
                            // Compute the index of face of mid-point flux and
                            // projection matrix.
                            const int idx_face_y = (j + 1) +
                                (i + 1)*(interior_dims[1] + 3);
                            
                            boost::multi_array<double, 2> W_array(
                                boost::extents[6][d_num_eqn]);
                            
                            /*
                             * Project primitive variables onto characteristic fields.
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
                                
                                std::vector<const double*> V_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_ptr.push_back(&Z_rho[si][idx_cell]);
                                }
                                V_ptr.push_back(&u[idx_cell]);
                                V_ptr.push_back(&v[idx_cell]);
                                V_ptr.push_back(&p[idx_cell]);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_ptr.push_back(&Z[si][idx_cell]);
                                }
                                
                                std::vector<double*> W_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_ptr.push_back(&W_array[m][ei]);
                                }
                                
                                projectPrimitiveVariablesToCharacteristicFields(
                                    W_ptr,
                                    V_ptr,
                                    R_y_inv_intercell);
                            }
                            
                            /*
                             * Do WENO interplation on characteristic variables to get W_B
                             * and W_T
                             */
                            
                            std::vector<double> W_B;
                            std::vector<double> W_T;
                            
                            performWENOInterpolation(W_B, W_T, W_array);
                            
                            /*
                             * Project characteristic variables back to physical fields.
                             */
                            
                            std::vector<double> Z_rho_B;
                            std::vector<double> Z_rho_T;
                            Z_rho_B.resize(d_num_species);
                            Z_rho_T.resize(d_num_species);
                            
                            std::vector<double> vel_B;
                            std::vector<double> vel_T;
                            vel_B.resize(d_dim.getValue());
                            vel_T.resize(d_dim.getValue());
                            
                            double p_B;
                            double p_T;
                            
                            std::vector<double> Z_B;
                            std::vector<double> Z_T;
                            Z_B.resize(d_num_species - 1);
                            Z_T.resize(d_num_species - 1);
                            
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
                            
                            std::vector<const double*> W_B_ptr;
                            std::vector<const double*> W_T_ptr;
                            
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                W_B_ptr.push_back(&W_B[ei]);
                                W_T_ptr.push_back(&W_T[ei]);
                            }
                            
                            std::vector<double*> V_B_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                V_B_ptr.push_back(&Z_rho_B[si]);
                            }
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                V_B_ptr.push_back(&vel_B[di]);
                            }
                            V_B_ptr.push_back(&p_B);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                V_B_ptr.push_back(&Z_B[si]);
                            }
                            
                            std::vector<double*> V_T_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                V_T_ptr.push_back(&Z_rho_T[si]);
                            }
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                V_T_ptr.push_back(&vel_T[di]);
                            }
                            V_T_ptr.push_back(&p_T);
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                V_T_ptr.push_back(&Z_T[si]);
                            }
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                V_B_ptr,
                                W_B_ptr,
                                R_y_intercell);
                            
                            projectCharacteristicVariablesToPhysicalFields(
                                V_T_ptr,
                                W_T_ptr,
                                R_y_intercell);
                            
                            /*
                             * Convert the primitive variables into conservative variables.
                             */
                            
                            std::vector<const double*> Z_rho_B_ptr;
                            std::vector<const double*> Z_rho_T_ptr;
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Z_rho_B_ptr.push_back(&Z_rho_B[si]);
                                Z_rho_T_ptr.push_back(&Z_rho_T[si]);
                            }
                            
                            double rho_B = d_equation_of_state->
                                getTotalDensity(
                                    Z_rho_B_ptr);
                            
                            double rho_T = d_equation_of_state->
                                getTotalDensity(
                                    Z_rho_T_ptr);
                            
                            std::vector<double> m_B;
                            std::vector<double> m_T;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_B.push_back(rho_B*vel_B[di]);
                                m_T.push_back(rho_T*vel_T[di]);
                            }
                            
                            std::vector<const double*> vel_B_ptr;
                            std::vector<const double*> vel_T_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                vel_B_ptr.push_back(&vel_B[di]);
                                vel_T_ptr.push_back(&vel_T[di]);
                            }
                            
                            std::vector<const double*> Z_B_ptr;
                            std::vector<const double*> Z_T_ptr;
                            for (int si = 0; si < d_num_species - 1; si++)
                            {
                                Z_B_ptr.push_back(&Z_B[si]);
                                Z_T_ptr.push_back(&Z_T[si]);
                            }
                            
                            double E_B = d_equation_of_state->
                                getTotalEnergyWithVolumeFraction(
                                    &rho_B,
                                    vel_B_ptr,
                                    &p_B,
                                    Z_B_ptr);
                            
                            double E_T = d_equation_of_state->
                                getTotalEnergyWithVolumeFraction(
                                    &rho_T,
                                    vel_T_ptr,
                                    &p_T,
                                    Z_T_ptr);
                            
                            /*
                             * If the WENO interpolated density or energy are negative,
                             * use constant interpolation.
                             */
                            if ((rho_B < 0) || (rho_T < 0) || (E_B < 0) || (E_T < 0))
                            {
                                // Compute the indices of bottom cell and top cell.
                                const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                const int idx_cell_T = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0];
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_B[si] = Z_rho[si][idx_cell_B];
                                    Z_rho_T[si] = Z_rho[si][idx_cell_T];
                                }
                                
                                m_B[0] = rho_u[idx_cell_B];
                                m_B[1] = rho_v[idx_cell_B];
                                m_T[0] = rho_u[idx_cell_T];
                                m_T[1] = rho_v[idx_cell_T];
                                
                                E_B = E[idx_cell_B];
                                E_T = E[idx_cell_T];
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_B[si] = Z[si][idx_cell_B];
                                    Z_T[si] = Z[si][idx_cell_T];
                                }
                            }
                            
                            /*
                             * Apply the Riemann solver.
                             */
                            
                            std::vector<const double*> m_B_ptr;
                            std::vector<const double*> m_T_ptr;
                            for (int di = 0; di < d_dim.getValue(); di++)
                            {
                                m_B_ptr.push_back(&m_B[di]);
                                m_T_ptr.push_back(&m_T[di]);
                            }
                            
                            std::vector<double*> F_y_midpoint_ptr;
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                F_y_midpoint_ptr.push_back(&F_y_midpoint[ei][idx_face_y]);
                            }
                            
                            std::vector<double*> vel_y_ptr;
                            for (int vi = 0; vi < d_dim.getValue(); vi++)
                            {
                                vel_y_ptr.push_back(&(velocity_intercell->getPointer(1, vi)[idx_face_y]));
                            }
                            
                            d_Riemann_solver.computeIntercellFluxAndVelocityForFiveEqnAllaire(
                                F_y_midpoint_ptr,
                                vel_y_ptr,
                                Z_rho_B_ptr,
                                Z_rho_T_ptr,
                                m_B_ptr,
                                m_T_ptr,
                                &E_B,
                                &E_T,
                                Z_B_ptr,
                                Z_T_ptr,
                                Y_DIRECTION);
                        }
                    }
                    
                    // Compute the flux in the x direction.
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0] + 1; i++)
                        {
                            // Compute the indices.
                            const int idx_face_x = i +
                                j*(interior_dims[0] + 1);
                            
                            const int idx_midpoint_x = (i + 1) +
                                (j + 1)*(interior_dims[0] + 3);
                            
                            const int idx_node_L = (i - 1 + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_node_R = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            // Compute the flux.
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                convective_flux->getPointer(0, ei)[idx_face_x] = dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x + 1] +
                                    F_x_midpoint[ei][idx_midpoint_x - 1]) -
                                    3.0/10*(F_x_node[ei][idx_node_R] +
                                    F_x_node[ei][idx_node_L]) +
                                    23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
                            }
                        }
                    }
                    
                    // Compute the flux in the y direction.
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int j = 0; j < interior_dims[1] + 1; j++)
                        {
                            // Compute the indices.
                            const int idx_face_y = j +
                                i*(interior_dims[1] + 1);
                            
                            const int idx_midpoint_y = (j + 1) +
                                (i + 1)*(interior_dims[1] + 3);
                            
                            const int idx_node_B = (i + d_num_ghosts[0]) +
                                (j - 1 + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            const int idx_node_T = (i + d_num_ghosts[0]) +
                                (j + d_num_ghosts[1])*ghostcell_dims[0];
                            
                            // Compute the flux.
                            for (int ei = 0; ei < d_num_eqn; ei++)
                            {
                                convective_flux->getPointer(1, ei)[idx_face_y] =
                                    dt*(1.0/30*(F_y_midpoint[ei][idx_midpoint_y + 1] +
                                    F_y_midpoint[ei][idx_midpoint_y - 1]) -
                                    3.0/10*(F_y_node[ei][idx_node_T] +
                                    F_y_node[ei][idx_node_B]) +
                                    23.0/15*F_y_midpoint[ei][idx_midpoint_y]);
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
                                
                                const int idx_face_x_L = (i + 1) +
                                    (j + 1)*(interior_dims[0] + 3);
                                
                                const int idx_face_x_R = (i + 2) +
                                    (j + 1)*(interior_dims[0] + 3);
                                
                                const int idx_face_y_B = (j + 1) +
                                    (i + 1)*(interior_dims[1] + 3);
                                
                                const int idx_face_y_T = (j + 2) +
                                    (i + 1)*(interior_dims[1] + 3);
                                
                                const int& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                                const int& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                                
                                const int& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                                const int& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                                
                                S[idx_cell_nghost] += dt*Z[si][idx_cell_wghost]*((u_R - u_L)/dx[0] + (v_T - v_B)/dx[1]);
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
                    std::vector<double*> F_x_node;
                    std::vector<double*> F_y_node;
                    std::vector<double*> F_z_node;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_node.push_back(convective_flux_node[0]->getPointer(ei));
                        F_y_node.push_back(convective_flux_node[1]->getPointer(ei));
                        F_z_node.push_back(convective_flux_node[2]->getPointer(ei));
                    }
                    std::vector<double*> F_x_midpoint;
                    std::vector<double*> F_y_midpoint;
                    std::vector<double*> F_z_midpoint;
                    for (int ei = 0; ei < d_num_eqn; ei++)
                    {
                        F_x_midpoint.push_back(convective_flux_midpoint->getPointer(0, ei));
                        F_y_midpoint.push_back(convective_flux_midpoint->getPointer(1, ei));
                        F_z_midpoint.push_back(convective_flux_midpoint->getPointer(2, ei));
                    }
                    
                    // Compute the field of total density, velocities, pressure and fluxes.
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
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_x_node[si][idx] = u[idx]*Z_rho[si][idx];
                                }
                                F_x_node[d_num_species][idx] = rho_u[idx]*u[idx] + p[idx];
                                F_x_node[d_num_species + 1][idx] = rho_u[idx]*v[idx];
                                F_x_node[d_num_species + 2][idx] = rho_u[idx]*w[idx];
                                F_x_node[d_num_species + 3][idx] = u[idx]*(E[idx] + p[idx]);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_x_node[d_num_species + 4 + si][idx] = u[idx]*Z[si][idx];
                                }
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_y_node[si][idx] = v[idx]*Z_rho[si][idx];
                                }
                                F_y_node[d_num_species][idx] = rho_v[idx]*u[idx];
                                F_y_node[d_num_species + 1][idx] = rho_v[idx]*v[idx] + p[idx];
                                F_y_node[d_num_species + 2][idx] = rho_v[idx]*w[idx];
                                F_y_node[d_num_species + 3][idx] = v[idx]*(E[idx] + p[idx]);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_y_node[d_num_species + 4 + si][idx] = v[idx]*Z[si][idx];
                                }
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    F_z_node[si][idx] = w[idx]*Z_rho[si][idx];
                                }
                                F_z_node[d_num_species][idx] = rho_w[idx]*u[idx];
                                F_z_node[d_num_species + 1][idx] = rho_w[idx]*v[idx];
                                F_z_node[d_num_species + 2][idx] = rho_w[idx]*w[idx] + p[idx];
                                F_z_node[d_num_species + 3][idx] = w[idx]*(E[idx] + p[idx]);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    F_z_node[d_num_species + 4 + si][idx] = w[idx]*Z[si][idx];
                                }
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
                            for (int i = -1; i < interior_dims[0] + 2; i++)
                            {
                                // Compute the indices of left cell, right cell and face of
                                // projection matrix.
                                const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_cell_R = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_face_x = (i + 1) +
                                    (j + 1)*(interior_dims[0] + 3) +
                                    (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                
                                // Get left and right quantities.
                                const double& rho_L = rho[idx_cell_L];
                                const double& rho_R = rho[idx_cell_R];
                                
                                const double& c_L = c[idx_cell_L];
                                const double& c_R = c[idx_cell_R];
                                
                                // Compute simply-averaged quantities.
                                const double rho_average = 0.5*(rho_L + rho_R);
                                const double c_average = 0.5*(c_L + c_R);
                                
                                std::vector<double> Z_rho_average;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_average.push_back(0.5*(Z_rho[si][idx_cell_L] +
                                        Z_rho[si][idx_cell_R]));
                                }
                                
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
                                            
                                        *R_x_intercell[ei][ej] = 0.0;
                                        *R_x_inv_intercell[ei][ej] = 0.0;
                                    }
                                }
                                
                                /*
                                 * Compute R_x_intercell.
                                 */
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    *R_x_intercell[si][0] = -0.5*Z_rho_average[si]/c_average;
                                    *R_x_intercell[si][si + 1] = 1.0;
                                    *R_x_intercell[si][d_num_eqn - 1] = 0.5*Z_rho_average[si]/c_average;
                                }
                                
                                *R_x_intercell[d_num_species][0] = 0.5;
                                *R_x_intercell[d_num_species][d_num_eqn - 1] = 0.5;
                                
                                *R_x_intercell[d_num_species + 1][d_num_species + 1] = 1.0;
                                
                                *R_x_intercell[d_num_species + 2][d_num_species + 2] = 1.0;
                                
                                *R_x_intercell[d_num_species + 3][0] = -0.5*rho_average*c_average;
                                *R_x_intercell[d_num_species + 3][d_num_eqn - 1] = 0.5*rho_average*c_average;
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    *R_x_intercell[d_num_species + 4 + si][d_num_species + 3 + si] = 1.0;
                                }
                                
                                /*
                                 * Compute R_x_inv_intercell.
                                 */
                                
                                *R_x_inv_intercell[0][d_num_species] = 1.0;
                                *R_x_inv_intercell[0][d_num_species + 3] = -1.0/(rho_average*c_average);
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    *R_x_inv_intercell[1 + si][si] = 1.0;
                                    *R_x_inv_intercell[1 + si][d_num_species + 3] = -Z_rho_average[si]/
                                        (rho_average*c_average*c_average);
                                }
                                
                                *R_x_inv_intercell[d_num_species + 1][d_num_species + 1] = 1.0;
                                
                                *R_x_inv_intercell[d_num_species + 2][d_num_species + 2] = 1.0;
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    *R_x_inv_intercell[d_num_species + 3 + si][d_num_species + 4 + si] = 1.0;
                                }
                                
                                *R_x_inv_intercell[2*d_num_species + 2][d_num_species] = 1.0;
                                *R_x_inv_intercell[2*d_num_species + 2][d_num_species + 3] = 1.0/(rho_average*c_average);
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
                            for (int j = -1; j < interior_dims[1] + 2; j++)
                            {
                                // Compute the indices of bottom cell, top cell and face of
                                // projection matrix.
                                const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_cell_T = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_face_y = (j + 1) +
                                    (k + 1)*(interior_dims[1] + 3) +
                                    (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                
                                // Get bottom and top quantities.
                                
                                const double& rho_B = rho[idx_cell_B];
                                const double& rho_T = rho[idx_cell_T];
                                
                                const double& c_B = c[idx_cell_B];
                                const double& c_T = c[idx_cell_T];
                                
                                // Compute simply-averaged quantities.
                                const double rho_average = 0.5*(rho_B + rho_T);
                                const double c_average = 0.5*(c_B + c_T);
                                
                                std::vector<double> Z_rho_average;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_average.push_back(0.5*(Z_rho[si][idx_cell_B] +
                                        Z_rho[si][idx_cell_T]));
                                }
                                
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
                                        
                                        *R_y_intercell[ei][ej] = 0.0;
                                        *R_y_inv_intercell[ei][ej] = 0.0;
                                    }
                                }
                                
                                /*
                                 * Compute R_y_intercell.
                                 */
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    *R_y_intercell[si][0] = -0.5*Z_rho_average[si]/c_average;
                                    *R_y_intercell[si][si + 1] = 1.0;
                                    *R_y_intercell[si][d_num_eqn - 1] = 0.5*Z_rho_average[si]/c_average;
                                }
                                
                                *R_y_intercell[d_num_species][d_num_species + 1] = 1.0;
                                
                                *R_y_intercell[d_num_species + 1][0] = 0.5;
                                *R_y_intercell[d_num_species + 1][d_num_eqn - 1] = 0.5;
                                
                                *R_y_intercell[d_num_species + 2][d_num_species + 2] = 1.0;
                                
                                *R_y_intercell[d_num_species + 3][0] = -0.5*rho_average*c_average;
                                *R_y_intercell[d_num_species + 3][d_num_eqn - 1] = 0.5*rho_average*c_average;
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    *R_y_intercell[d_num_species + 4 + si][d_num_species + 3 + si] = 1.0;
                                }
                                
                                /*
                                 * Compute R_y_inv_intercell.
                                 */
                                
                                *R_y_inv_intercell[0][d_num_species + 1] = 1.0;
                                *R_y_inv_intercell[0][d_num_species + 3] = -1.0/(rho_average*c_average);
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    *R_y_inv_intercell[1 + si][si] = 1.0;
                                    *R_y_inv_intercell[1 + si][d_num_species + 3] = -Z_rho_average[si]/
                                        (rho_average*c_average*c_average);
                                }
                                
                                *R_y_inv_intercell[d_num_species + 1][d_num_species] = 1.0;
                                
                                *R_y_inv_intercell[d_num_species + 2][d_num_species + 2] = 1.0;
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    *R_y_inv_intercell[d_num_species + 3 + si][d_num_species + 4 + si] = 1.0;
                                }
                                
                                *R_y_inv_intercell[2*d_num_species + 2][d_num_species + 1] = 1.0;
                                *R_y_inv_intercell[2*d_num_species + 2][d_num_species + 3] = 1.0/(rho_average*c_average);
                            }
                        }
                    }
                    
                    /*
                     * Compute the projection matrix and its inverse at the face normal to the
                     * z direction.
                     */
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int k = -1; k < interior_dims[2] + 2; k++)
                            {
                                // Compute the indices of back cell, front cell and face of
                                // projection matrix.
                                const int idx_cell_B = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_cell_F = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_face_z = (k + 1) +
                                    (i + 1)*(interior_dims[2] + 3) +
                                    (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                
                                // Get back and front quantities.
                                
                                const double& rho_B = rho[idx_cell_B];
                                const double& rho_F = rho[idx_cell_F];
                                
                                const double& c_B = c[idx_cell_B];
                                const double& c_F = c[idx_cell_F];
                                
                                // Compute simply-averaged quantities.
                                const double rho_average = 0.5*(rho_B + rho_F);
                                const double c_average = 0.5*(c_B + c_F);
                                
                                std::vector<double> Z_rho_average;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_average.push_back(0.5*(Z_rho[si][idx_cell_B] +
                                        Z_rho[si][idx_cell_F]));
                                }
                                
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
                                        
                                        *R_z_intercell[ei][ej] = 0.0;
                                        *R_z_inv_intercell[ei][ej] = 0.0;
                                    }
                                }
                                
                                /*
                                 * Compute R_z_intercell.
                                 */
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    *R_z_intercell[si][0] = -0.5*Z_rho_average[si]/c_average;
                                    *R_z_intercell[si][si + 1] = 1.0;
                                    *R_z_intercell[si][d_num_eqn - 1] = 0.5*Z_rho_average[si]/c_average;
                                }
                                
                                *R_z_intercell[d_num_species][d_num_species + 1] = 1.0;
                                
                                *R_z_intercell[d_num_species + 1][d_num_species + 2] = 1.0;
                                
                                *R_z_intercell[d_num_species + 2][0] = 0.5;
                                *R_z_intercell[d_num_species + 2][d_num_eqn - 1] = 0.5;
                                
                                *R_z_intercell[d_num_species + 3][0] = -0.5*rho_average*c_average;
                                *R_z_intercell[d_num_species + 3][d_num_eqn - 1] = 0.5*rho_average*c_average;
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    *R_z_intercell[d_num_species + 4 + si][d_num_species + 3 + si] = 1.0;
                                }
                                
                                /*
                                 * Compute R_z_inv_intercell.
                                 */
                                
                                *R_z_inv_intercell[0][d_num_species + 2] = 1.0;
                                *R_z_inv_intercell[0][d_num_species + 3] = -1.0/(rho_average*c_average);
                                
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    *R_z_inv_intercell[1 + si][si] = 1.0;
                                    *R_z_inv_intercell[1 + si][d_num_species + 3] = -Z_rho_average[si]/
                                        (rho_average*c_average*c_average);
                                }
                                
                                *R_z_inv_intercell[d_num_species + 1][d_num_species] = 1.0;
                                
                                *R_z_inv_intercell[d_num_species + 2][d_num_species + 1] = 1.0;
                                
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    *R_z_inv_intercell[d_num_species + 3 + si][d_num_species + 4 + si] = 1.0;
                                }
                                
                                *R_z_inv_intercell[2*d_num_species + 2][d_num_species + 2] = 1.0;
                                *R_z_inv_intercell[2*d_num_species + 2][d_num_species + 3] = 1.0/(rho_average*c_average);
                            }
                        }
                    }
                    
                    // Compute the mid-point flux in the x direction.
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = -1; i < interior_dims[0] + 2; i++)
                            {
                                // Compute the index of face of mid-point flux and
                                // projection matrix.
                                const int idx_face_x = (i + 1) +
                                    (j + 1)*(interior_dims[0] + 3) +
                                    (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                
                                boost::multi_array<double, 2> W_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                /*
                                 * Project primitive variables onto characteristic fields.
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
                                    
                                    std::vector<const double*> V_ptr;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ptr.push_back(&Z_rho[si][idx_cell]);
                                    }
                                    V_ptr.push_back(&u[idx_cell]);
                                    V_ptr.push_back(&v[idx_cell]);
                                    V_ptr.push_back(&w[idx_cell]);
                                    V_ptr.push_back(&p[idx_cell]);
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        V_ptr.push_back(&Z[si][idx_cell]);
                                    }
                                    
                                    std::vector<double*> W_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        W_ptr.push_back(&W_array[m][ei]);
                                    }
                                    projectPrimitiveVariablesToCharacteristicFields(
                                        W_ptr,
                                        V_ptr,
                                        R_x_inv_intercell);
                                }
                                
                                /*
                                 * Do WENO interplation on characteristic variables to get W_L
                                 * and W_R
                                 */
                                
                                std::vector<double> W_L;
                                std::vector<double> W_R;
                                
                                performWENOInterpolation(W_L, W_R, W_array);
                                
                                /*
                                 * Project characteristic variables back to physical fields.
                                 */
                                
                                std::vector<double> Z_rho_L;
                                std::vector<double> Z_rho_R;
                                Z_rho_L.resize(d_num_species);
                                Z_rho_R.resize(d_num_species);
                                
                                std::vector<double> vel_L;
                                std::vector<double> vel_R;
                                vel_L.resize(d_dim.getValue());
                                vel_R.resize(d_dim.getValue());
                                
                                double p_L;
                                double p_R;
                                
                                std::vector<double> Z_L;
                                std::vector<double> Z_R;
                                Z_L.resize(d_num_species - 1);
                                Z_R.resize(d_num_species - 1);
                                
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
                                
                                std::vector<const double*> W_L_ptr;
                                std::vector<const double*> W_R_ptr;
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_L_ptr.push_back(&W_L[ei]);
                                    W_R_ptr.push_back(&W_R[ei]);
                                }
                                
                                std::vector<double*> V_L_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_L_ptr.push_back(&Z_rho_L[si]);
                                }
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_L_ptr.push_back(&vel_L[di]);
                                }
                                V_L_ptr.push_back(&p_L);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_L_ptr.push_back(&Z_L[si]);
                                }
                                
                                std::vector<double*> V_R_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_R_ptr.push_back(&Z_rho_R[si]);
                                }
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_R_ptr.push_back(&vel_R[di]);
                                }
                                V_R_ptr.push_back(&p_R);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_R_ptr.push_back(&Z_R[si]);
                                }
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_L_ptr,
                                    W_L_ptr,
                                    R_x_intercell);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_R_ptr,
                                    W_R_ptr,
                                    R_x_intercell);
                                
                                /*
                                 * Convert the primitive variables into conservative variables.
                                 */
                                
                                std::vector<const double*> Z_rho_L_ptr;
                                std::vector<const double*> Z_rho_R_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_L_ptr.push_back(&Z_rho_L[si]);
                                    Z_rho_R_ptr.push_back(&Z_rho_R[si]);
                                }
                                
                                double rho_L = d_equation_of_state->
                                    getTotalDensity(
                                        Z_rho_L_ptr);
                                
                                double rho_R = d_equation_of_state->
                                    getTotalDensity(
                                        Z_rho_R_ptr);
                                
                                std::vector<double> m_L;
                                std::vector<double> m_R;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_L.push_back(rho_L*vel_L[di]);
                                    m_R.push_back(rho_R*vel_R[di]);
                                }
                                
                                std::vector<const double*> vel_L_ptr;
                                std::vector<const double*> vel_R_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    vel_L_ptr.push_back(&vel_L[di]);
                                    vel_R_ptr.push_back(&vel_R[di]);
                                }
                                
                                std::vector<const double*> Z_L_ptr;
                                std::vector<const double*> Z_R_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_L_ptr.push_back(&Z_L[si]);
                                    Z_R_ptr.push_back(&Z_R[si]);
                                }
                                
                                double E_L = d_equation_of_state->
                                    getTotalEnergyWithVolumeFraction(
                                        &rho_L,
                                        vel_L_ptr,
                                        &p_L,
                                        Z_L_ptr);
                                
                                double E_R = d_equation_of_state->
                                    getTotalEnergyWithVolumeFraction(
                                        &rho_R,
                                        vel_R_ptr,
                                        &p_R,
                                        Z_R_ptr);
                                
                                /*
                                 * If the WENO interpolated density or energy are negative,
                                 * use constant interpolation.
                                 */
                                if ((rho_L < 0) || (rho_R < 0) || (E_L < 0) || (E_R < 0))
                                {
                                    // Compute the indices of left cell and right cell.
                                    const int idx_cell_L = (i - 1 + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_R = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho_L[si] = Z_rho[si][idx_cell_L];
                                        Z_rho_R[si] = Z_rho[si][idx_cell_R];
                                    }
                                    
                                    m_L[0] = rho_u[idx_cell_L];
                                    m_L[1] = rho_v[idx_cell_L];
                                    m_L[2] = rho_w[idx_cell_L];
                                    m_R[0] = rho_u[idx_cell_R];
                                    m_R[1] = rho_v[idx_cell_R];
                                    m_R[2] = rho_w[idx_cell_R];
                                    
                                    E_L = E[idx_cell_L];
                                    E_R = E[idx_cell_R];
                                    
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        Z_L[si] = Z[si][idx_cell_L];
                                        Z_R[si] = Z[si][idx_cell_R];
                                    }
                                }
                                
                                /*
                                 * Apply the Riemann solver.
                                 */
                                
                                std::vector<const double*> m_L_ptr;
                                std::vector<const double*> m_R_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_L_ptr.push_back(&m_L[di]);
                                    m_R_ptr.push_back(&m_R[di]);
                                }
                                
                                std::vector<double*> F_x_midpoint_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_x_midpoint_ptr.push_back(&F_x_midpoint[ei][idx_face_x]);
                                }
                                
                                std::vector<double*> vel_x_ptr;
                                for (int vi = 0; vi < d_dim.getValue(); vi++)
                                {
                                    vel_x_ptr.push_back(&(velocity_intercell->getPointer(0, vi)[idx_face_x]));
                                }
                                
                                d_Riemann_solver.computeIntercellFluxAndVelocityForFiveEqnAllaire(
                                    F_x_midpoint_ptr,
                                    vel_x_ptr,
                                    Z_rho_L_ptr,
                                    Z_rho_R_ptr,
                                    m_L_ptr,
                                    m_R_ptr,
                                    &E_L,
                                    &E_R,
                                    Z_L_ptr,
                                    Z_R_ptr,
                                    X_DIRECTION);
                            }
                        }
                    }
                    
                    // Compute the mid-point flux in the y direction.
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = -1; j < interior_dims[1] + 2; j++)
                            {
                                // Compute the index of face of mid-point flux and
                                // projection matrix.
                                const int idx_face_y = (j + 1) +
                                    (k + 1)*(interior_dims[1] + 3) +
                                    (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                
                                boost::multi_array<double, 2> W_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                /*
                                 * Project primitive variables onto characteristic fields.
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
                                    
                                    std::vector<const double*> V_ptr;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ptr.push_back(&Z_rho[si][idx_cell]);
                                    }
                                    V_ptr.push_back(&u[idx_cell]);
                                    V_ptr.push_back(&v[idx_cell]);
                                    V_ptr.push_back(&w[idx_cell]);
                                    V_ptr.push_back(&p[idx_cell]);
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        V_ptr.push_back(&Z[si][idx_cell]);
                                    }
                                    
                                    std::vector<double*> W_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        W_ptr.push_back(&W_array[m][ei]);
                                    }
                                    
                                    projectPrimitiveVariablesToCharacteristicFields(
                                        W_ptr,
                                        V_ptr,
                                        R_y_inv_intercell);
                                }
                                
                                /*
                                 * Do WENO interplation on characteristic variables to get W_B
                                 * and W_T
                                 */
                                
                                std::vector<double> W_B;
                                std::vector<double> W_T;
                                
                                performWENOInterpolation(W_B, W_T, W_array);
                                
                                /*
                                 * Project characteristic variables back to physical fields.
                                 */
                                
                                std::vector<double> Z_rho_B;
                                std::vector<double> Z_rho_T;
                                Z_rho_B.resize(d_num_species);
                                Z_rho_T.resize(d_num_species);
                                
                                std::vector<double> vel_B;
                                std::vector<double> vel_T;
                                vel_B.resize(d_dim.getValue());
                                vel_T.resize(d_dim.getValue());
                                
                                double p_B;
                                double p_T;
                                
                                std::vector<double> Z_B;
                                std::vector<double> Z_T;
                                Z_B.resize(d_num_species - 1);
                                Z_T.resize(d_num_species - 1);
                                
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
                                
                                std::vector<const double*> W_B_ptr;
                                std::vector<const double*> W_T_ptr;
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_B_ptr.push_back(&W_B[ei]);
                                    W_T_ptr.push_back(&W_T[ei]);
                                }
                                
                                std::vector<double*> V_B_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_B_ptr.push_back(&Z_rho_B[si]);
                                }
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_B_ptr.push_back(&vel_B[di]);
                                }
                                V_B_ptr.push_back(&p_B);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_B_ptr.push_back(&Z_B[si]);
                                }
                                
                                std::vector<double*> V_T_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_T_ptr.push_back(&Z_rho_T[si]);
                                }
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_T_ptr.push_back(&vel_T[di]);
                                }
                                V_T_ptr.push_back(&p_T);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_T_ptr.push_back(&Z_T[si]);
                                }
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_B_ptr,
                                    W_B_ptr,
                                    R_y_intercell);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_T_ptr,
                                    W_T_ptr,
                                    R_y_intercell);
                                
                                /*
                                 * Convert the primitive variables into conservative variables.
                                 */
                                
                                std::vector<const double*> Z_rho_B_ptr;
                                std::vector<const double*> Z_rho_T_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_B_ptr.push_back(&Z_rho_B[si]);
                                    Z_rho_T_ptr.push_back(&Z_rho_T[si]);
                                }
                                
                                double rho_B = d_equation_of_state->
                                    getTotalDensity(
                                        Z_rho_B_ptr);
                                
                                double rho_T = d_equation_of_state->
                                    getTotalDensity(
                                        Z_rho_T_ptr);
                                
                                std::vector<double> m_B;
                                std::vector<double> m_T;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_B.push_back(rho_B*vel_B[di]);
                                    m_T.push_back(rho_T*vel_T[di]);
                                }
                                
                                std::vector<const double*> vel_B_ptr;
                                std::vector<const double*> vel_T_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    vel_B_ptr.push_back(&vel_B[di]);
                                    vel_T_ptr.push_back(&vel_T[di]);
                                }
                                
                                std::vector<const double*> Z_B_ptr;
                                std::vector<const double*> Z_T_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_B_ptr.push_back(&Z_B[si]);
                                    Z_T_ptr.push_back(&Z_T[si]);
                                }
                                
                                double E_B = d_equation_of_state->
                                    getTotalEnergyWithVolumeFraction(
                                        &rho_B,
                                        vel_B_ptr,
                                        &p_B,
                                        Z_B_ptr);
                                
                                double E_T = d_equation_of_state->
                                    getTotalEnergyWithVolumeFraction(
                                        &rho_T,
                                        vel_T_ptr,
                                        &p_T,
                                        Z_T_ptr);
                                
                                /*
                                 * If the WENO interpolated density or energy are negative,
                                 * use constant interpolation.
                                 */
                                if ((rho_B < 0) || (rho_T < 0) || (E_B < 0) || (E_T < 0))
                                {
                                    // Compute the indices of bottom cell and top cell.
                                    const int idx_cell_B = (i + d_num_ghosts[0]) +
                                        (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_T = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho_B[si] = Z_rho[si][idx_cell_B];
                                        Z_rho_T[si] = Z_rho[si][idx_cell_T];
                                    }
                                    
                                    m_B[0] = rho_u[idx_cell_B];
                                    m_B[1] = rho_v[idx_cell_B];
                                    m_B[2] = rho_w[idx_cell_B];
                                    m_T[0] = rho_u[idx_cell_T];
                                    m_T[1] = rho_v[idx_cell_T];
                                    m_T[2] = rho_w[idx_cell_T];
                                    
                                    E_B = E[idx_cell_B];
                                    E_T = E[idx_cell_T];
                                    
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        Z_B[si] = Z[si][idx_cell_B];
                                        Z_T[si] = Z[si][idx_cell_T];
                                    }
                                }
                                
                                /*
                                 * Apply the Riemann solver.
                                 */
                                
                                std::vector<const double*> m_B_ptr;
                                std::vector<const double*> m_T_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_B_ptr.push_back(&m_B[di]);
                                    m_T_ptr.push_back(&m_T[di]);
                                }
                                
                                std::vector<double*> F_y_midpoint_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_y_midpoint_ptr.push_back(&F_y_midpoint[ei][idx_face_y]);
                                }
                                
                                std::vector<double*> vel_y_ptr;
                                for (int vi = 0; vi < d_dim.getValue(); vi++)
                                {
                                    vel_y_ptr.push_back(&(velocity_intercell->getPointer(1, vi)[idx_face_y]));
                                }
                                
                                d_Riemann_solver.computeIntercellFluxAndVelocityForFiveEqnAllaire(
                                    F_y_midpoint_ptr,
                                    vel_y_ptr,
                                    Z_rho_B_ptr,
                                    Z_rho_T_ptr,
                                    m_B_ptr,
                                    m_T_ptr,
                                    &E_B,
                                    &E_T,
                                    Z_B_ptr,
                                    Z_T_ptr,
                                    Y_DIRECTION);
                            }
                        }
                    }
                    
                    // Compute the mid-point flux in the z direction.
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int k = -1; k < interior_dims[2] + 2; k++)
                            {
                                // Compute the index of face of mid-point flux and
                                // projection matrix.
                                const int idx_face_z = (k + 1) +
                                    (i + 1)*(interior_dims[2] + 3) +
                                    (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                
                                boost::multi_array<double, 2> W_array(
                                    boost::extents[6][d_num_eqn]);
                                
                                /*
                                 * Project primitive variables onto characteristic fields.
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
                                    
                                    std::vector<const double*> V_ptr;
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        V_ptr.push_back(&Z_rho[si][idx_cell]);
                                    }
                                    V_ptr.push_back(&u[idx_cell]);
                                    V_ptr.push_back(&v[idx_cell]);
                                    V_ptr.push_back(&w[idx_cell]);
                                    V_ptr.push_back(&p[idx_cell]);
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        V_ptr.push_back(&Z[si][idx_cell]);
                                    }
                                    
                                    std::vector<double*> W_ptr;
                                    for (int ei = 0; ei < d_num_eqn; ei++)
                                    {
                                        W_ptr.push_back(&W_array[m][ei]);
                                    }
                                    
                                    projectPrimitiveVariablesToCharacteristicFields(
                                        W_ptr,
                                        V_ptr,
                                        R_z_inv_intercell);
                                }
                                
                                /*
                                 * Do WENO interplation on characteristic variables to get W_B
                                 * and W_F
                                 */
                                
                                std::vector<double> W_B;
                                std::vector<double> W_F;
                                
                                performWENOInterpolation(W_B, W_F, W_array);
                                
                                /*
                                 * Project characteristic variables back to physical fields.
                                 */
                                
                                std::vector<double> Z_rho_B;
                                std::vector<double> Z_rho_F;
                                Z_rho_B.resize(d_num_species);
                                Z_rho_F.resize(d_num_species);
                                
                                std::vector<double> vel_B;
                                std::vector<double> vel_F;
                                vel_B.resize(d_dim.getValue());
                                vel_F.resize(d_dim.getValue());
                                
                                double p_B;
                                double p_F;
                                
                                std::vector<double> Z_B;
                                std::vector<double> Z_F;
                                Z_B.resize(d_num_species - 1);
                                Z_F.resize(d_num_species - 1);
                                
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
                                
                                std::vector<const double*> W_B_ptr;
                                std::vector<const double*> W_F_ptr;
                                
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    W_B_ptr.push_back(&W_B[ei]);
                                    W_F_ptr.push_back(&W_F[ei]);
                                }
                                
                                std::vector<double*> V_B_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_B_ptr.push_back(&Z_rho_B[si]);
                                }
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_B_ptr.push_back(&vel_B[di]);
                                }
                                V_B_ptr.push_back(&p_B);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_B_ptr.push_back(&Z_B[si]);
                                }
                                
                                std::vector<double*> V_F_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    V_F_ptr.push_back(&Z_rho_F[si]);
                                }
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    V_F_ptr.push_back(&vel_F[di]);
                                }
                                V_F_ptr.push_back(&p_F);
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    V_F_ptr.push_back(&Z_F[si]);
                                }
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_B_ptr,
                                    W_B_ptr,
                                    R_z_intercell);
                                
                                projectCharacteristicVariablesToPhysicalFields(
                                    V_F_ptr,
                                    W_F_ptr,
                                    R_z_intercell);
                                
                                /*
                                 * Convert the primitive variables into conservative variables.
                                 */
                                
                                std::vector<const double*> Z_rho_B_ptr;
                                std::vector<const double*> Z_rho_F_ptr;
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Z_rho_B_ptr.push_back(&Z_rho_B[si]);
                                    Z_rho_F_ptr.push_back(&Z_rho_F[si]);
                                }
                                
                                double rho_B = d_equation_of_state->
                                    getTotalDensity(
                                        Z_rho_B_ptr);
                                
                                double rho_F = d_equation_of_state->
                                    getTotalDensity(
                                        Z_rho_F_ptr);
                                
                                std::vector<double> m_B;
                                std::vector<double> m_F;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_B.push_back(rho_B*vel_B[di]);
                                    m_F.push_back(rho_F*vel_F[di]);
                                }
                                
                                std::vector<const double*> vel_B_ptr;
                                std::vector<const double*> vel_F_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    vel_B_ptr.push_back(&vel_B[di]);
                                    vel_F_ptr.push_back(&vel_F[di]);
                                }
                                
                                std::vector<const double*> Z_B_ptr;
                                std::vector<const double*> Z_F_ptr;
                                for (int si = 0; si < d_num_species - 1; si++)
                                {
                                    Z_B_ptr.push_back(&Z_B[si]);
                                    Z_F_ptr.push_back(&Z_F[si]);
                                }
                                
                                double E_B = d_equation_of_state->
                                    getTotalEnergyWithVolumeFraction(
                                        &rho_B,
                                        vel_B_ptr,
                                        &p_B,
                                        Z_B_ptr);
                                
                                double E_F = d_equation_of_state->
                                    getTotalEnergyWithVolumeFraction(
                                        &rho_F,
                                        vel_F_ptr,
                                        &p_F,
                                        Z_F_ptr);
                                
                                /*
                                 * If the WENO interpolated density or energy are negative,
                                 * use constant interpolation.
                                 */
                                if ((rho_B < 0) || (rho_F < 0) || (E_B < 0) || (E_F < 0))
                                {
                                    // Compute the indices of back cell and front cell.
                                    const int idx_cell_B = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                    
                                    const int idx_cell_F = (i + d_num_ghosts[0]) +
                                        (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                        (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                        
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Z_rho_B[si] = Z_rho[si][idx_cell_B];
                                        Z_rho_F[si] = Z_rho[si][idx_cell_F];
                                    }
                                    
                                    m_B[0] = rho_u[idx_cell_B];
                                    m_B[1] = rho_v[idx_cell_B];
                                    m_B[2] = rho_w[idx_cell_B];
                                    m_F[0] = rho_u[idx_cell_F];
                                    m_F[1] = rho_v[idx_cell_F];
                                    m_F[2] = rho_w[idx_cell_F];
                                    
                                    E_B = E[idx_cell_B];
                                    E_F = E[idx_cell_F];
                                    
                                    for (int si = 0; si < d_num_species - 1; si++)
                                    {
                                        Z_B[si] = Z[si][idx_cell_B];
                                        Z_F[si] = Z[si][idx_cell_F];
                                    }
                                }
                                
                                /*
                                 * Apply the Riemann solver.
                                 */
                                
                                std::vector<const double*> m_B_ptr;
                                std::vector<const double*> m_F_ptr;
                                for (int di = 0; di < d_dim.getValue(); di++)
                                {
                                    m_B_ptr.push_back(&m_B[di]);
                                    m_F_ptr.push_back(&m_F[di]);
                                }
                                
                                std::vector<double*> F_z_midpoint_ptr;
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    F_z_midpoint_ptr.push_back(&F_z_midpoint[ei][idx_face_z]);
                                }
                                
                                std::vector<double*> vel_z_ptr;
                                for (int vi = 0; vi < d_dim.getValue(); vi++)
                                {
                                    vel_z_ptr.push_back(&(velocity_intercell->getPointer(2, vi)[idx_face_z]));
                                }
                                
                                d_Riemann_solver.computeIntercellFluxAndVelocityForFiveEqnAllaire(
                                    F_z_midpoint_ptr,
                                    vel_z_ptr,
                                    Z_rho_B_ptr,
                                    Z_rho_F_ptr,
                                    m_B_ptr,
                                    m_F_ptr,
                                    &E_B,
                                    &E_F,
                                    Z_B_ptr,
                                    Z_F_ptr,
                                    Z_DIRECTION);
                            }
                        }
                    }
                    
                    // Compute the flux in the x direction.
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0] + 1; i++)
                            {
                                // Compute the indices.
                                const int idx_face_x = i +
                                    j*(interior_dims[0] + 1) +
                                    k*(interior_dims[0] + 1)*interior_dims[1];
                                
                                const int idx_midpoint_x = (i + 1) +
                                    (j + 1)*(interior_dims[0] + 3) +
                                    (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                
                                const int idx_node_L = (i - 1 + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_node_R = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                // Compute the flux.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    convective_flux->getPointer(0, ei)[idx_face_x] = dt*(1.0/30*(F_x_midpoint[ei][idx_midpoint_x + 1] +
                                        F_x_midpoint[ei][idx_midpoint_x - 1]) -
                                        3.0/10*(F_x_node[ei][idx_node_R] +
                                        F_x_node[ei][idx_node_L]) +
                                        23.0/15*F_x_midpoint[ei][idx_midpoint_x]);
                                }
                            }
                        }
                    }
                    
                    // Compute the flux in the y direction.
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        for (int k = 0; k < interior_dims[2]; k++)
                        {
                            for (int j = 0; j < interior_dims[1] + 1; j++)
                            {
                                // Compute the indices.
                                const int idx_face_y = j +
                                    k*(interior_dims[1] + 1) +
                                    i*(interior_dims[1] + 1)*interior_dims[2];
                                
                                const int idx_midpoint_y = (j + 1) +
                                    (k + 1)*(interior_dims[1] + 3) +
                                    (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                
                                const int idx_node_B = (i + d_num_ghosts[0]) +
                                    (j - 1 + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_node_T = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                // Compute the flux.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    convective_flux->getPointer(1, ei)[idx_face_y] =
                                        dt*(1.0/30*(F_y_midpoint[ei][idx_midpoint_y + 1] +
                                        F_y_midpoint[ei][idx_midpoint_y - 1]) -
                                        3.0/10*(F_y_node[ei][idx_node_T] +
                                        F_y_node[ei][idx_node_B]) +
                                        23.0/15*F_y_midpoint[ei][idx_midpoint_y]);
                                }
                            }
                        }
                    }
                    
                    // Compute the flux in the z direction.
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            for (int k = 0; k < interior_dims[2] + 1; k++)
                            {
                                // Compute the indices.
                                const int idx_face_z = k +
                                    i*(interior_dims[2] + 1) +
                                    j*(interior_dims[2] + 1)*interior_dims[0];
                                
                                const int idx_midpoint_z = (k + 1) +
                                    (i + 1)*(interior_dims[2] + 3) +
                                    (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                
                                const int idx_node_B = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k - 1 + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                const int idx_node_F = (i + d_num_ghosts[0]) +
                                    (j + d_num_ghosts[1])*ghostcell_dims[0] +
                                    (k + d_num_ghosts[2])*ghostcell_dims[0]*ghostcell_dims[1];
                                
                                // Compute the flux.
                                for (int ei = 0; ei < d_num_eqn; ei++)
                                {
                                    convective_flux->getPointer(2, ei)[idx_face_z] =
                                        dt*(1.0/30*(F_z_midpoint[ei][idx_midpoint_z + 1] +
                                        F_z_midpoint[ei][idx_midpoint_z - 1]) -
                                        3.0/10*(F_z_node[ei][idx_node_F] +
                                        F_z_node[ei][idx_node_B]) +
                                        23.0/15*F_z_midpoint[ei][idx_midpoint_z]);
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
                                    
                                    const int idx_face_x_L = (i + 1) +
                                        (j + 1)*(interior_dims[0] + 3) +
                                        (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                    
                                    const int idx_face_x_R = (i + 2) +
                                        (j + 1)*(interior_dims[0] + 3) +
                                        (k + 1)*(interior_dims[0] + 3)*(interior_dims[1] + 2);
                                    
                                    const int idx_face_y_B = (j + 1) +
                                        (k + 1)*(interior_dims[1] + 3) +
                                        (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                    
                                    const int idx_face_y_T = (j + 2) +
                                        (k + 1)*(interior_dims[1] + 3) +
                                        (i + 1)*(interior_dims[1] + 3)*(interior_dims[2] + 2);
                                    
                                    const int idx_face_z_B = (k + 1) +
                                        (i + 1)*(interior_dims[2] + 3) +
                                        (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                    
                                    const int idx_face_z_F = (k + 2) +
                                        (i + 1)*(interior_dims[2] + 3) +
                                        (j + 1)*(interior_dims[2] + 3)*(interior_dims[0] + 2);
                                    
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
 * Convert primitive variables into characteristic variables.
 */
void
ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL::projectPrimitiveVariablesToCharacteristicFields(
    std::vector<double*> characteristic_variables,
    const std::vector<const double*> primitive_variables,
    const boost::multi_array<const double*, 2> projection_matrix)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(characteristic_variables.size() ==
        primitive_variables.size());
    
    TBOX_ASSERT(projection_matrix.num_dimensions() == 2);
    
    TBOX_ASSERT(projection_matrix.shape()[0] ==
        projection_matrix.shape()[1]);
    
    TBOX_ASSERT(projection_matrix.size() ==
        primitive_variables.size());
#endif
    
    std::vector<double*>&                         W               = characteristic_variables;
    const std::vector<const double*>&             V               = primitive_variables;
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
 * Convert characteristic variables into primitive variables.
 */
void
ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL::projectCharacteristicVariablesToPhysicalFields(
    std::vector<double*> primitive_variables,
    const std::vector<const double*> characteristic_variables,
    const boost::multi_array<const double*, 2> projection_matrix_inv)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(characteristic_variables.size() ==
        primitive_variables.size());
    
    TBOX_ASSERT(projection_matrix_inv.num_dimensions() == 2);
    
    TBOX_ASSERT(projection_matrix_inv.shape()[0] ==
        projection_matrix_inv.shape()[1]);
    
    TBOX_ASSERT(projection_matrix_inv.size() ==
        primitive_variables.size());
#endif
    
    std::vector<double*>&                         V           = primitive_variables;
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
 * Compute beta's.
 */
boost::multi_array<double, 2>
ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL::computeBeta(
    const boost::multi_array<double, 2>& W_array)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
    TBOX_ASSERT(static_cast<int>(W_array.shape()[1]) == d_num_eqn);
#endif
    
    boost::multi_array<double, 2> beta(
        boost::extents[d_num_eqn][4]);
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        beta[ei][0] = 1.0/3*(W_array[0][ei]*(4*W_array[0][ei] - 19*W_array[1][ei] + 11*W_array[2][ei]) +
            W_array[1][ei]*(25*W_array[1][ei] - 31*W_array[2][ei]) + 10*W_array[2][ei]*W_array[2][ei]);
        
        beta[ei][1] = 1.0/3*(W_array[1][ei]*(4*W_array[1][ei] - 13*W_array[2][ei] + 5*W_array[3][ei]) +
            13*W_array[2][ei]*(W_array[2][ei] - W_array[3][ei]) + 4*W_array[3][ei]*W_array[3][ei]);
        
        beta[ei][2] = 1.0/3*(W_array[2][ei]*(10*W_array[2][ei] - 31*W_array[3][ei] + 11*W_array[4][ei]) +
            W_array[3][ei]*(25*W_array[3][ei] - 19*W_array[4][ei]) + 4*W_array[4][ei]*W_array[4][ei]);
    }
    
    return beta;
}


/*
 * Compute beta_tilde's.
 */
boost::multi_array<double, 2>
ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL::computeBetaTilde(
    const boost::multi_array<double, 2>& W_array)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(static_cast<int>(W_array.shape()[0]) == 6);
    TBOX_ASSERT(static_cast<int>(W_array.shape()[1]) == d_num_eqn);
#endif
    
    boost::multi_array<double, 2> beta_tilde(
        boost::extents[d_num_eqn][4]);
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        beta_tilde[ei][0] = 1.0/3*(W_array[5][ei]*(4*W_array[5][ei] - 19*W_array[4][ei] + 11*W_array[3][ei]) +
            W_array[4][ei]*(25*W_array[4][ei] - 31*W_array[3][ei]) + 10*W_array[3][ei]*W_array[3][ei]);
        
        beta_tilde[ei][1] = 1.0/3*(W_array[4][ei]*(4*W_array[4][ei] - 13*W_array[3][ei] + 5*W_array[2][ei]) +
            13*W_array[3][ei]*(W_array[3][ei] - W_array[2][ei]) + 4*W_array[2][ei]*W_array[2][ei]);
        
        beta_tilde[ei][2] = 1.0/3*(W_array[3][ei]*(10*W_array[3][ei] - 31*W_array[2][ei] + 11*W_array[1][ei]) +
            W_array[2][ei]*(25*W_array[2][ei] - 19*W_array[1][ei]) + 4*W_array[1][ei]*W_array[1][ei]);
    }
    
    return beta_tilde;
}


/*
 * Perform WENO interpolation.
 */
void
ConvectiveFluxReconstructorWCNS_JS5_HLLC_HLL::performWENOInterpolation(
    std::vector<double>& W_L,
    std::vector<double>& W_R,
    const boost::multi_array<double, 2>& W_array)
{
    /*
     * Compute the beta's.
     */
    
    boost::multi_array<double, 2> beta = computeBeta(W_array);
    boost::multi_array<double, 2> beta_tilde = computeBetaTilde(W_array);
    
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        /*
         * Do WENO interpolation on current characteristic variables.
         */
        
        std::vector<double> W_L_r;
        std::vector<double> W_R_r;
        
        // Do the interpolation from different stencils.
        for (int r = 0; r < 3; r++)
        {
            W_L_r.push_back(0.0);
            W_R_r.push_back(0.0);
            
            for (int m = r; m < 3 + r; m++)
            {
                W_L_r[r] += d_weights_c[r][m - r]*W_array[m][ei];
                W_R_r[r] += d_weights_c[r][m - r]*W_array[6 - m - 1][ei];
            }
        }
        
        /*
         * Compute W_L of the current characteristic variable.
         */
        
        // Compute the weights alpha.
        std::vector<double> alpha;
        for (int r = 0; r < 3; r++)
        {
            alpha.push_back(d_weights_d[r]/pow((beta[ei][r] + EPSILON), 2));
        }
        
        // Sum up the weights alpha.
        double alpha_sum = 0.0;
        for (int r = 0; r < 3; r++)
        {
            alpha_sum += alpha[r];
        }
        
        // Compute the nonlinear weights omega.
        std::vector<double> omega;
        for (int r = 0; r < 3; r++)
        {
            omega.push_back(alpha[r]/alpha_sum);
        }
        
        // Compute the W_L.
        W_L.push_back(0.0);
        for (int r = 0; r < 3; r++)
        {
            W_L[ei] += omega[r]*W_L_r[r];
        }
        
        /*
         * Compute W_R of the current characteristic variable.
         */
        
        // Compute the weights alpha_tilde.
        std::vector<double> alpha_tilde;
        for (int r = 0; r < 3; r++)
        {
            alpha_tilde.push_back(d_weights_d[r]/pow((beta_tilde[ei][r] + EPSILON), 2));
        }
        
        // Sum up the weights alpha_tilde.
        double alpha_tilde_sum = 0.0;
        for (int r = 0; r < 3; r++)
        {
            alpha_tilde_sum += alpha_tilde[r];
        }
        
        // Compute the nonlinear weights omega_tilde.
        std::vector<double> omega_tilde;
        for (int r = 0; r < 3; r++)
        {
            omega_tilde.push_back(alpha_tilde[r]/alpha_tilde_sum);
        }
        
        // Compute the W_R.
        W_R.push_back(0.0);
        for (int r = 0; r < 3; r++)
        {
            W_R[ei] += omega_tilde[r]*W_R_r[r];
        }
    }
}
