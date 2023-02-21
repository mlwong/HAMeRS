#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <map>

NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const HAMERS_SHARED_PTR<tbox::Database>& nonconservative_diffusive_flux_divergence_operator_db):
        NonconservativeDiffusiveFluxDivergenceOperator(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model,
            nonconservative_diffusive_flux_divergence_operator_db)
{
    d_num_diff_ghosts = hier::IntVector::getOne(d_dim)*6;
}


/*
 * Print all characteristics of the non-conservative diffusive flux divergence operator class.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::printClassData(
    std::ostream& os) const
{
    os << "\nPrint NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder object..."
       << std::endl;
    
    os << std::endl;
    
    os << "NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder: this = "
       << (NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Put the characteristics of the non-conservative diffusive flux divergence operator class into
 * the restart database.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::putToRestart(
   const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    restart_db->putString("d_nonconservative_diffusive_flux_divergence_operator", "TWELFTH_ORDER");
}


/*
 * Compute the first derivatives in the x-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::computeFirstDerivativesInX(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_x,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_x_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_x,
    const std::vector<std::vector<int> >& data_component_idx_x,
    const hier::Patch& patch)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_x[ei].size()) ==
                    static_cast<int>(data_component_idx_x[ei].size()));
    }
#endif
    
    derivative_x.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dx_0_inv = double(1)/dx[0];
    
    const double a_n =  double(6)/double(7);
    const double b_n = -double(15)/double(56);
    const double c_n =  double(5)/double(63);
    const double d_n = -double(1)/double(56);
    const double e_n =  double(3)/double(1155);
    const double f_n = -double(1)/double(5544);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudx = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + num_diff_ghosts_0;
                        
                        const int idx_data_LLLLLL = i - 6 + num_subghosts_0_data;
                        const int idx_data_LLLLL  = i - 5 + num_subghosts_0_data;
                        const int idx_data_LLLL   = i - 4 + num_subghosts_0_data;
                        const int idx_data_LLL    = i - 3 + num_subghosts_0_data;
                        const int idx_data_LL     = i - 2 + num_subghosts_0_data;
                        const int idx_data_L      = i - 1 + num_subghosts_0_data;
                        const int idx_data_R      = i + 1 + num_subghosts_0_data;
                        const int idx_data_RR     = i + 2 + num_subghosts_0_data;
                        const int idx_data_RRR    = i + 3 + num_subghosts_0_data;
                        const int idx_data_RRRR   = i + 4 + num_subghosts_0_data;
                        const int idx_data_RRRRR  = i + 5 + num_subghosts_0_data;
                        const int idx_data_RRRRRR = i + 6 + num_subghosts_0_data;
                        
                        dudx[idx] = (a_n*(u[idx_data_R] - u[idx_data_L]) +
                                     b_n*(u[idx_data_RR] - u[idx_data_LL]) +
                                     c_n*(u[idx_data_RRR] - u[idx_data_LLL]) +
                                     d_n*(u[idx_data_RRRR] - u[idx_data_LLLL]) +
                                     e_n*(u[idx_data_RRRRR] - u[idx_data_LLLLL]) +
                                     f_n*(u[idx_data_RRRRRR] - u[idx_data_LLLLLL]))*
                                        dx_0_inv;
                    }
                    
                    std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudx = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_x[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    for (int j = -3; j < interior_dim_1 + 3; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_LLLLLL = (i - 6 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LLLLL = (i - 5 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LLLL = (i - 4 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRRR = (i + 4 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRRRR = (i + 5 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRRRRR = (i + 6 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            dudx[idx] = (a_n*(u[idx_data_R] - u[idx_data_L]) +
                                         b_n*(u[idx_data_RR] - u[idx_data_LL]) +
                                         c_n*(u[idx_data_RRR] - u[idx_data_LLL]) +
                                         d_n*(u[idx_data_RRRR] - u[idx_data_LLLL]) +
                                         e_n*(u[idx_data_RRRRR] - u[idx_data_LLLLL]) +
                                         f_n*(u[idx_data_RRRRRR] - u[idx_data_LLLLLL]))*
                                            dx_0_inv;
                        }
                    }
                    
                    std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudx = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_x[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = -3; k < interior_dim_2 + 3; k++)
                    {
                        for (int j = -3; j < interior_dim_1 + 3; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_LLLLLL = (i - 6 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LLLLL = (i - 5 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LLLL = (i - 4 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRRR = (i + 4 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRRRR = (i + 5 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRRRRR = (i + 6 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                dudx[idx] = (a_n*(u[idx_data_R] - u[idx_data_L]) +
                                             b_n*(u[idx_data_RR] - u[idx_data_LL]) +
                                             c_n*(u[idx_data_RRR] - u[idx_data_LLL]) +
                                             d_n*(u[idx_data_RRRR] - u[idx_data_LLLL]) +
                                             e_n*(u[idx_data_RRRRR] - u[idx_data_LLLLL]) +
                                             f_n*(u[idx_data_RRRRRR] - u[idx_data_LLLLLL]))*
                                                dx_0_inv;
                            }
                        }
                    }
                    
                    std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the first derivatives in the y-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::computeFirstDerivativesInY(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_y,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_y_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_y,
    const std::vector<std::vector<int> >& data_component_idx_y,
    const hier::Patch& patch)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_y[ei].size()) ==
                    static_cast<int>(data_component_idx_y[ei].size()));
    }
#endif
    
    derivative_y.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dx_1_inv = double(1)/dx[1];
    
    const double a_n =  double(6)/double(7);
    const double b_n = -double(15)/double(56);
    const double c_n =  double(5)/double(63);
    const double d_n = -double(1)/double(56);
    const double e_n =  double(3)/double(1155);
    const double f_n = -double(1)/double(5544);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::"
            << "computeFirstDerivativesInY()\n"
            << "There isn't y-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            derivative_y[ei].reserve(static_cast<int>(data_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_y[ei][vi];
                
                if (derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                    == derivative_y_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_y[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudy = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_y[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_y[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = -3; i < interior_dim_0 + 3; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_BBBBBB = (i + num_subghosts_0_data) +
                                (j - 6 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BBBBB = (i + num_subghosts_0_data) +
                                (j - 5 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BBBB = (i + num_subghosts_0_data) +
                                (j - 4 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BBB = (i + num_subghosts_0_data) +
                                (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BB = (i + num_subghosts_0_data) +
                                (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_B = (i + num_subghosts_0_data) +
                                (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_T = (i + num_subghosts_0_data) +
                                (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TT = (i + num_subghosts_0_data) +
                                (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTT = (i + num_subghosts_0_data) +
                                (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTTT = (i + num_subghosts_0_data) +
                                (j + 4 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTTTT = (i + num_subghosts_0_data) +
                                (j + 5 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTTTTT = (i + num_subghosts_0_data) +
                                (j + 6 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            dudy[idx] = (a_n*(u[idx_data_T] - u[idx_data_B]) +
                                         b_n*(u[idx_data_TT] - u[idx_data_BB]) +
                                         c_n*(u[idx_data_TTT] - u[idx_data_BBB]) +
                                         d_n*(u[idx_data_TTTT] - u[idx_data_BBBB]) +
                                         e_n*(u[idx_data_TTTTT] - u[idx_data_BBBBB]) +
                                         f_n*(u[idx_data_TTTTTT] - u[idx_data_BBBBBB]))*
                                            dx_1_inv;
                        }
                    }
                    
                    std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_y_computed.insert(derivative_pair);
                }
                
                derivative_y[ei].push_back(
                    derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            derivative_y[ei].reserve(static_cast<int>(data_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_y[ei][vi];
                
                if (derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                    == derivative_y_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_y[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudy = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_y[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_y[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = -3; k < interior_dim_2 + 3; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = -3; i < interior_dim_0 + 3; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBBBBB = (i + num_subghosts_0_data) +
                                    (j - 6 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BBBBB = (i + num_subghosts_0_data) +
                                    (j - 5 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BBBB = (i + num_subghosts_0_data) +
                                    (j - 4 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_T = (i + num_subghosts_0_data) +
                                    (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TT = (i + num_subghosts_0_data) +
                                    (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTT = (i + num_subghosts_0_data) +
                                    (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTTT = (i + num_subghosts_0_data) +
                                    (j + 4 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTTTT = (i + num_subghosts_0_data) +
                                    (j + 5 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTTTTT = (i + num_subghosts_0_data) +
                                    (j + 6 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                dudy[idx] = (a_n*(u[idx_data_T] - u[idx_data_B]) +
                                             b_n*(u[idx_data_TT] - u[idx_data_BB]) +
                                             c_n*(u[idx_data_TTT] - u[idx_data_BBB]) +
                                             d_n*(u[idx_data_TTTT] - u[idx_data_BBBB]) +
                                             e_n*(u[idx_data_TTTTT] - u[idx_data_BBBBB]) +
                                             f_n*(u[idx_data_TTTTTT] - u[idx_data_BBBBBB]))*
                                                dx_1_inv;
                            }
                        }
                    }
                    
                    std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_y_computed.insert(derivative_pair);
                }
                
                derivative_y[ei].push_back(
                    derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the first derivatives in the z-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::computeFirstDerivativesInZ(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_z,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_z_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_z,
    const std::vector<std::vector<int> >& data_component_idx_z,
    const hier::Patch& patch)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_z[ei].size()) ==
                    static_cast<int>(data_component_idx_z[ei].size()));
    }
#endif
    
    derivative_z.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dx_2_inv = double(1)/dx[2];
    
    const double a_n =  double(6)/double(7);
    const double b_n = -double(15)/double(56);
    const double c_n =  double(5)/double(63);
    const double d_n = -double(1)/double(56);
    const double e_n =  double(3)/double(1155);
    const double f_n = -double(1)/double(5544);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::"
            << "computeFirstDerivativesInZ()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::"
            << "computeFirstDerivativesInZ()\n"
            << "There isn't z-direction for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            derivative_z[ei].reserve(static_cast<int>(data_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_z[ei][vi];
                
                if (derivative_z_computed.find(data_z[ei][vi]->getPointer(u_idx))
                    == derivative_z_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_z[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudz = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_z[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_z[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = -3; j < interior_dim_1 + 3; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = -3; i < interior_dim_0 + 3; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBBBBB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 6 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BBBBB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 5 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BBBB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 4 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_F = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FFF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FFFF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 4 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FFFFF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 5 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FFFFFF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 6 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                dudz[idx] = (a_n*(u[idx_data_F] - u[idx_data_B]) +
                                             b_n*(u[idx_data_FF] - u[idx_data_BB]) +
                                             c_n*(u[idx_data_FFF] - u[idx_data_BBB]) +
                                             d_n*(u[idx_data_FFFF] - u[idx_data_BBBB]) +
                                             e_n*(u[idx_data_FFFFF] - u[idx_data_BBBBB]) +
                                             f_n*(u[idx_data_FFFFFF] - u[idx_data_BBBBBB]))*
                                                dx_2_inv;
                            }
                        }
                    }
                    
                    std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_z_computed.insert(derivative_pair);
                }
                
                derivative_z[ei].push_back(
                    derivative_z_computed.find(data_z[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the second derivatives in the x-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::computeSecondDerivativesInX(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_x,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_x_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_x,
    const std::vector<std::vector<int> >& data_component_idx_x,
    const hier::Patch& patch)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_x[ei].size()) ==
                    static_cast<int>(data_component_idx_x[ei].size()));
    }
#endif
    
    derivative_x.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dx_sq_inv = double(1)/(dx[0]*dx[0]);
    
    const double a_n = -double(5369)/double(1800);
    const double b_n =  double(12)/double(7);
    const double c_n = -double(15)/double(56);
    const double d_n =  double(10)/double(189);
    const double e_n = -double(1)/double(112);
    const double f_n =  double(2)/double(1925);
    const double g_n = -double(1)/double(16632);
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udx2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + num_diff_ghosts_0;
                        
                        const int idx_data_LLLLLL = i - 6 + num_subghosts_0_data;
                        const int idx_data_LLLLL  = i - 5 + num_subghosts_0_data;
                        const int idx_data_LLLL   = i - 4 + num_subghosts_0_data;
                        const int idx_data_LLL    = i - 3 + num_subghosts_0_data;
                        const int idx_data_LL     = i - 2 + num_subghosts_0_data;
                        const int idx_data_L      = i - 1 + num_subghosts_0_data;
                        const int idx_data        = i     + num_subghosts_0_data;
                        const int idx_data_R      = i + 1 + num_subghosts_0_data;
                        const int idx_data_RR     = i + 2 + num_subghosts_0_data;
                        const int idx_data_RRR    = i + 3 + num_subghosts_0_data;
                        const int idx_data_RRRR   = i + 4 + num_subghosts_0_data;
                        const int idx_data_RRRRR  = i + 5 + num_subghosts_0_data;
                        const int idx_data_RRRRRR = i + 6 + num_subghosts_0_data;
                        
                        d2udx2[idx] = (a_n*u[idx_data] +
                                       b_n*(u[idx_data_L] + u[idx_data_R]) +
                                       c_n*(u[idx_data_LL] + u[idx_data_RR]) +
                                       d_n*(u[idx_data_LLL] + u[idx_data_RRR]) +
                                       e_n*(u[idx_data_LLLL] + u[idx_data_RRRR]) +
                                       f_n*(u[idx_data_LLLLL] + u[idx_data_RRRRR]) +
                                       g_n*(u[idx_data_LLLLLL] + u[idx_data_RRRRRR]))*
                                        dx_sq_inv;
                    }
                    
                    std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udx2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_x[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_LLLLLL = (i - 6 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LLLLL = (i - 5 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LLLL = (i - 4 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRRR = (i + 4 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRRRR = (i + 5 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRRRRR = (i + 6 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            d2udx2[idx] = (a_n*u[idx_data] +
                                           b_n*(u[idx_data_L] + u[idx_data_R]) +
                                           c_n*(u[idx_data_LL] + u[idx_data_RR]) +
                                           d_n*(u[idx_data_LLL] + u[idx_data_RRR]) +
                                           e_n*(u[idx_data_LLLL] + u[idx_data_RRRR]) +
                                           f_n*(u[idx_data_LLLLL] + u[idx_data_RRRRR]) +
                                           g_n*(u[idx_data_LLLLLL] + u[idx_data_RRRRRR]))*
                                            dx_sq_inv;
                        }
                    }
                    
                    std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udx2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_x[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_LLLLLL = (i - 6 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LLLLL = (i - 5 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LLLL = (i - 4 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRRR = (i + 4 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRRRR = (i + 5 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRRRRR = (i + 6 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                d2udx2[idx] = (a_n*u[idx_data] +
                                               b_n*(u[idx_data_L] + u[idx_data_R]) +
                                               c_n*(u[idx_data_LL] + u[idx_data_RR]) +
                                               d_n*(u[idx_data_LLL] + u[idx_data_RRR]) +
                                               e_n*(u[idx_data_LLLL] + u[idx_data_RRRR]) +
                                               f_n*(u[idx_data_LLLLL] + u[idx_data_RRRRR]) +
                                               g_n*(u[idx_data_LLLLLL] + u[idx_data_RRRRRR]))*
                                                dx_sq_inv;
                            }
                        }
                    }
                    
                    std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the second derivatives in the y-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::computeSecondDerivativesInY(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_y,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_y_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_y,
    const std::vector<std::vector<int> >& data_component_idx_y,
    const hier::Patch& patch)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_y[ei].size()) ==
                    static_cast<int>(data_component_idx_y[ei].size()));
    }
#endif
    
    derivative_y.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dy_sq_inv = double(1)/(dx[1]*dx[1]);
    
    const double a_n = -double(5369)/double(1800);
    const double b_n =  double(12)/double(7);
    const double c_n = -double(15)/double(56);
    const double d_n =  double(10)/double(189);
    const double e_n = -double(1)/double(112);
    const double f_n =  double(2)/double(1925);
    const double g_n = -double(1)/double(16632);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::"
            << "computeSecondDerivativesInY()\n"
            << "There isn't y-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            derivative_y[ei].reserve(static_cast<int>(data_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_y[ei][vi];
                
                if (derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                    == derivative_y_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_y[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udy2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_y[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_y[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_BBBBBB = (i + num_subghosts_0_data) +
                                (j - 6 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BBBBB = (i + num_subghosts_0_data) +
                                (j - 5 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BBBB = (i + num_subghosts_0_data) +
                                (j - 4 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BBB = (i + num_subghosts_0_data) +
                                (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BB = (i + num_subghosts_0_data) +
                                (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_B = (i + num_subghosts_0_data) +
                                (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_T = (i + num_subghosts_0_data) +
                                (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TT = (i + num_subghosts_0_data) +
                                (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTT = (i + num_subghosts_0_data) +
                                (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTTT = (i + num_subghosts_0_data) +
                                (j + 4 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTTTT = (i + num_subghosts_0_data) +
                                (j + 5 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTTTTT = (i + num_subghosts_0_data) +
                                (j + 6 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            d2udy2[idx] = (a_n*u[idx_data] +
                                           b_n*(u[idx_data_B] + u[idx_data_T]) +
                                           c_n*(u[idx_data_BB] + u[idx_data_TT]) +
                                           d_n*(u[idx_data_BBB] + u[idx_data_TTT]) +
                                           e_n*(u[idx_data_BBBB] + u[idx_data_TTTT]) +
                                           f_n*(u[idx_data_BBBBB] + u[idx_data_TTTTT]) +
                                           g_n*(u[idx_data_BBBBBB] + u[idx_data_TTTTTT]))*
                                            dy_sq_inv;
                        }
                    }
                    
                    std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_y_computed.insert(derivative_pair);
                }
                
                derivative_y[ei].push_back(
                    derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            derivative_y[ei].reserve(static_cast<int>(data_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_y[ei][vi];
                
                if (derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                    == derivative_y_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_y[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udy2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_y[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_y[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBBBBB = (i + num_subghosts_0_data) +
                                    (j - 6 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BBBBB = (i + num_subghosts_0_data) +
                                    (j - 5 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BBBB = (i + num_subghosts_0_data) +
                                    (j - 4 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_T = (i + num_subghosts_0_data) +
                                    (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TT = (i + num_subghosts_0_data) +
                                    (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTT = (i + num_subghosts_0_data) +
                                    (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTTT = (i + num_subghosts_0_data) +
                                    (j + 4 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTTTT = (i + num_subghosts_0_data) +
                                    (j + 5 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTTTTT = (i + num_subghosts_0_data) +
                                    (j + 6 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                d2udy2[idx] = (a_n*u[idx_data] +
                                               b_n*(u[idx_data_B] + u[idx_data_T]) +
                                               c_n*(u[idx_data_BB] + u[idx_data_TT]) +
                                               d_n*(u[idx_data_BBB] + u[idx_data_TTT]) +
                                               e_n*(u[idx_data_BBBB] + u[idx_data_TTTT]) +
                                               f_n*(u[idx_data_BBBBB] + u[idx_data_TTTTT]) +
                                               g_n*(u[idx_data_BBBBBB] + u[idx_data_TTTTTT]))*
                                                dy_sq_inv;
                            }
                        }
                    }
                    
                    std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_y_computed.insert(derivative_pair);
                }
                
                derivative_y[ei].push_back(
                    derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the second derivatives in the z-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::computeSecondDerivativesInZ(
    std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& derivative_z,
    std::map<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > >& derivative_z_computed,
    const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > > >& data_z,
    const std::vector<std::vector<int> >& data_component_idx_z,
    const hier::Patch& patch)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        TBOX_ASSERT(static_cast<int>(data_z[ei].size()) ==
                    static_cast<int>(data_component_idx_z[ei].size()));
    }
#endif
    
    derivative_z.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dz_sq_inv = double(1)/(dx[2]*dx[2]);
    
    const double a_n = -double(5369)/double(1800);
    const double b_n =  double(12)/double(7);
    const double c_n = -double(15)/double(56);
    const double d_n =  double(10)/double(189);
    const double e_n = -double(1)/double(112);
    const double f_n =  double(2)/double(1925);
    const double g_n = -double(1)/double(16632);
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::"
            << "computeSecondDerivativesInZ()\n"
            << "There isn't z-direction for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorTwelfthOrder::"
            << "computeSecondDerivativesInZ()\n"
            << "There isn't z-direction for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            derivative_z[ei].reserve(static_cast<int>(data_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_z[ei][vi];
                
                if (derivative_z_computed.find(data_z[ei][vi]->getPointer(u_idx))
                    == derivative_z_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_z[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    HAMERS_SHARED_PTR<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udz2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_z[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_z[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBBBBB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 6 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BBBBB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 5 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BBBB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 4 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_F = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FFF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FFFF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 4 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FFFFF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 5 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FFFFFF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 6 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                d2udz2[idx] = (a_n*u[idx_data] +
                                               b_n*(u[idx_data_B] + u[idx_data_F]) +
                                               c_n*(u[idx_data_BB] + u[idx_data_FF]) +
                                               d_n*(u[idx_data_BBB] + u[idx_data_FFF]) +
                                               e_n*(u[idx_data_BBBB] + u[idx_data_FFFF]) +
                                               f_n*(u[idx_data_BBBBB] + u[idx_data_FFFFF]) +
                                               g_n*(u[idx_data_BBBBBB] + u[idx_data_FFFFFF]))*
                                                dz_sq_inv;
                            }
                        }
                    }
                    
                    std::pair<double*, HAMERS_SHARED_PTR<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_z_computed.insert(derivative_pair);
                }
                
                derivative_z[ei].push_back(
                    derivative_z_computed.find(data_z[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}
