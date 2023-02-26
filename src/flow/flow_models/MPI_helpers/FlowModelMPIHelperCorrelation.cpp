#include "flow/flow_models/MPI_helpers/FlowModelMPIHelperCorrelation.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"

/*
 * Compute correlation with only x-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<std::vector<Real> >& averaged_quantities,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<Real> correlation;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousXDirection():\n"
            << "Not implemented for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        std::vector<const Real*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        Real* corr_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        correlation.resize(finest_level_dim_0);
        Real* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i]  = Real(0);
            corr_global[i] = Real(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const Real weight = Real(dx[1])/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                Real corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((Real) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        d_mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        std::vector<const Real*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        Real* corr_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        correlation.resize(finest_level_dim_0);
        Real* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i]  = Real(0);
            corr_global[i] = Real(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const Real weight = Real(dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    Real corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((Real) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        d_mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/*
 * Compute correlation with only x-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousXDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const std::vector<std::vector<Real> >& averaged_quantities,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<Real> correlation;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousXDirection():\n"
            << "Not implemented for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        std::vector<const Real*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        Real* corr_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        correlation.resize(finest_level_dim_0);
        Real* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i]  = Real(0);
            corr_global[i] = Real(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const Real weight = Real(dx[1])/L_y;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                Real corr = Real(0);
                                if (use_reciprocal[0])
                                {
                                    corr = Real(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                }
                                else
                                {
                                    corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                }
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        corr *= (Real(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    else
                                    {
                                        corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((Real) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        d_mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = d_finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        std::vector<const Real*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        Real* corr_local = (Real*)std::malloc(finest_level_dim_0*sizeof(Real));
        
        correlation.resize(finest_level_dim_0);
        Real* corr_global = correlation.data();
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            corr_local[i]  = Real(0);
            corr_global[i] = Real(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_0 = ratio_to_finest_level[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const Real weight = Real(dx[1]*dx[2])/(L_y*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int ii = 0; ii < ratio_to_finest_level_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratio_to_finest_level_0 + ii;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    Real corr = Real(0);
                                    if (use_reciprocal[0])
                                    {
                                        corr = Real(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                    }
                                    else
                                    {
                                        corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    }
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            corr *= (Real(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                        }
                                        else
                                        {
                                            corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                        }
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((Real) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        d_mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_0,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/*
 * Compute correlation with only y-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousYDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<std::vector<Real> >& averaged_quantities,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<Real> correlation;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousYDirection():\n"
            << "Cannot compute averaged value with only y-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
        std::vector<const Real*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        Real* corr_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        correlation.resize(finest_level_dim_1);
        Real* corr_global = correlation.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            corr_local[j]  = Real(0);
            corr_global[j] = Real(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const Real weight = Real(dx[0])/L_x;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                Real corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((Real) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        d_mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_1,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        std::vector<const Real*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        Real* corr_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        correlation.resize(finest_level_dim_1);
        Real* corr_global = correlation.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            corr_local[j]  = Real(0);
            corr_global[j] = Real(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const Real weight = Real(dx[0]*dx[2])/(L_x*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    Real corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((Real) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        d_mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_1,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/*
 * Compute correlation with only y-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousYDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const std::vector<std::vector<Real> >& averaged_quantities,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<Real> correlation;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousYDirection():\n"
            << "Cannot compute averaged value with only y-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        
        std::vector<const Real*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        Real* corr_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        correlation.resize(finest_level_dim_1);
        Real* corr_global = correlation.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            corr_local[j]  = Real(0);
            corr_global[j] = Real(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                }
                
                const Real weight = Real(dx[0])/L_x;
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j);
                            
                            int n_overlapped = 1;
                            
                            for (hier::BoxContainer::BoxContainerConstIterator iob(
                                    patch_overlapped_visible_boxes.begin());
                                 iob != patch_overlapped_visible_boxes.end();
                                 iob++)
                            {
                                const hier::Box& patch_overlapped_visible_box = *iob;
                                
                                if (patch_overlapped_visible_box.contains(idx_pt))
                                {
                                    n_overlapped++;
                                }
                            }
                            
                            /*
                             * Compute the linear indices and the data to add.
                             */
                            
                            for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                
                                const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0];
                                
                                Real corr = Real(0);
                                if (use_reciprocal[0])
                                {
                                    corr = Real(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                }
                                else
                                {
                                    corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                }
                                
                                for (int qi = 1; qi < num_quantities; qi++)
                                {
                                    const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi];
                                    
                                    if (use_reciprocal[qi])
                                    {
                                        corr *= (Real(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    else
                                    {
                                        corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                }
                                
                                corr_local[idx_fine] += (corr*weight/((Real) n_overlapped));
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        d_mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_1,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(corr_local);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = d_finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_z = Real(x_hi[2] - x_lo[2]);
        
        std::vector<const Real*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        Real* corr_local = (Real*)std::malloc(finest_level_dim_1*sizeof(Real));
        
        correlation.resize(finest_level_dim_1);
        Real* corr_global = correlation.data();
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            corr_local[j]  = Real(0);
            corr_global[j] = Real(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_1 = ratio_to_finest_level[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const Real weight = Real(dx[0]*dx[2])/(L_x*L_z);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int jj = 0; jj < ratio_to_finest_level_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratio_to_finest_level_1 + jj;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    Real corr = Real(0);
                                    if (use_reciprocal[0])
                                    {
                                        corr = Real(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                    }
                                    else
                                    {
                                        corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    }
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            corr *= (Real(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                        }
                                        else
                                        {
                                            corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                        }
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((Real) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        d_mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_1,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/*
 * Compute correlation with only z-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousZDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<std::vector<Real> >& averaged_quantities,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<Real> correlation;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_2 = d_finest_level_dims[2];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        std::vector<const Real*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        Real* corr_local = (Real*)std::malloc(finest_level_dim_2*sizeof(Real));
        
        correlation.resize(finest_level_dim_2);
        Real* corr_global = correlation.data();
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            corr_local[k]  = Real(0);
            corr_global[k] = Real(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_2 = ratio_to_finest_level[2];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int kk = 0; kk < ratio_to_finest_level_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratio_to_finest_level_2 + kk;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    Real corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((Real) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        d_mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_2,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}


/*
 * Compute correlation with only z-direction as inhomogeneous direction.
 */
std::vector<Real>
FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousZDirection(
    const std::vector<std::string>& quantity_names,
    const std::vector<int>& component_indices,
    const std::vector<bool>& use_reciprocal,
    const std::vector<std::vector<Real> >& averaged_quantities,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    int num_quantities = static_cast<int>(quantity_names.size());
    
    TBOX_ASSERT(static_cast<int>(component_indices.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(use_reciprocal.size()) == num_quantities);
    TBOX_ASSERT(static_cast<int>(averaged_quantities.size()) == num_quantities);
    
    std::vector<Real> correlation;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *d_patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelMPIHelperCorrelation::getQuantityCorrelationWithInhomogeneousZDirection():\n"
            << "Cannot compute averaged value with only z-direction as inhomogeneous direction for"
            << " two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_2 = d_finest_level_dims[2];
        
        /*
         * Get the size of the physical domain.
         */
        
        const Real L_x = Real(x_hi[0] - x_lo[0]);
        const Real L_y = Real(x_hi[1] - x_lo[1]);
        
        std::vector<const Real*> u_qi_avg_global;
        u_qi_avg_global.reserve(num_quantities);
        for (int qi = 0; qi < num_quantities; qi++)
        {
            u_qi_avg_global.push_back(averaged_quantities[qi].data());
        }
        
        Real* corr_local = (Real*)std::malloc(finest_level_dim_2*sizeof(Real));
        
        correlation.resize(finest_level_dim_2);
        Real* corr_global = correlation.data();
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            corr_local[k]  = Real(0);
            corr_global[k] = Real(0);
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratio_to_coarsest_level =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratio_to_coarsest_level *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratio_to_finest_level = d_ratio_finest_level_to_coarsest_level/ratio_to_coarsest_level;
            
            const int ratio_to_finest_level_2 = ratio_to_finest_level[2];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and data in the flow model and compute the corresponding
                 * correlation.
                 */
                
                setupFlowModelAndRegisterPatchWithDataContext(*patch, data_context);
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    num_subghosts_of_data.insert(
                        std::pair<std::string, hier::IntVector>(quantity_names[qi], hier::IntVector::getZero(d_dim)));
                }
                
                registerDerivedVariables(num_subghosts_of_data);
                
                allocateMemoryForDerivedCellData();
                
                computeDerivedCellData();
                
                /*
                 * Get the pointers to data inside the flow model.
                 */
                
                std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > > data_quantities;
                data_quantities.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    data_quantities[qi] = getCellData(quantity_names[qi]);
                }
                
                std::vector<Real*> u_qi;
                u_qi.resize(num_quantities);
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    u_qi[qi] = data_quantities[qi]->getPointer(component_indices[qi]);
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                std::vector<int> num_ghosts_0_u_qi;
                std::vector<int> num_ghosts_1_u_qi;
                std::vector<int> num_ghosts_2_u_qi;
                std::vector<int> ghostcell_dim_0_u_qi;
                std::vector<int> ghostcell_dim_1_u_qi;
                num_ghosts_0_u_qi.reserve(num_quantities);
                num_ghosts_1_u_qi.reserve(num_quantities);
                num_ghosts_2_u_qi.reserve(num_quantities);
                ghostcell_dim_0_u_qi.reserve(num_quantities);
                ghostcell_dim_1_u_qi.reserve(num_quantities);
                
                for (int qi = 0; qi < num_quantities; qi++)
                {
                    const hier::IntVector num_ghosts_u_qi = data_quantities[qi]->getGhostCellWidth();
                    const hier::IntVector ghostcell_dims_u_qi = data_quantities[qi]->getGhostBox().numberCells();
                    
                    num_ghosts_0_u_qi.push_back(num_ghosts_u_qi[0]);
                    num_ghosts_1_u_qi.push_back(num_ghosts_u_qi[1]);
                    num_ghosts_2_u_qi.push_back(num_ghosts_u_qi[2]);
                    ghostcell_dim_0_u_qi.push_back(ghostcell_dims_u_qi[0]);
                    ghostcell_dim_1_u_qi.push_back(ghostcell_dims_u_qi[1]);
                }
                
                const Real weight = Real(dx[0]*dx[1])/(L_x*L_y);
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    const int interior_dim_1 = interior_dims[1];
                    const int interior_dim_2 = interior_dims[2];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int idx_lo_1 = index_lo[1];
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the index of the data point and count how many times the data is repeated.
                                 */
                                
                                const hier::Index idx_pt(idx_lo_0 + i, idx_lo_1 + j, idx_lo_2 + k);
                                
                                int n_overlapped = 1;
                                
                                for (hier::BoxContainer::BoxContainerConstIterator iob(
                                        patch_overlapped_visible_boxes.begin());
                                     iob != patch_overlapped_visible_boxes.end();
                                     iob++)
                                {
                                    const hier::Box& patch_overlapped_visible_box = *iob;
                                    
                                    if (patch_overlapped_visible_box.contains(idx_pt))
                                    {
                                        n_overlapped++;
                                    }
                                }
                                
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                for (int kk = 0; kk < ratio_to_finest_level_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratio_to_finest_level_2 + kk;
                                    
                                    const int idx_q0 = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[0]) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_u_qi[0])*ghostcell_dim_0_u_qi[0] +
                                        (relative_idx_lo_2 + k + num_ghosts_2_u_qi[0])*ghostcell_dim_0_u_qi[0]*
                                            ghostcell_dim_1_u_qi[0];
                                    
                                    Real corr = Real(0);
                                    if (use_reciprocal[0])
                                    {
                                        corr = Real(1)/(u_qi[0][idx_q0]) - u_qi_avg_global[0][idx_fine];
                                    }
                                    else
                                    {
                                        corr = u_qi[0][idx_q0] - u_qi_avg_global[0][idx_fine];
                                    }
                                    
                                    for (int qi = 1; qi < num_quantities; qi++)
                                    {
                                        const int idx_qi = (relative_idx_lo_0 + i + num_ghosts_0_u_qi[qi]) +
                                            (relative_idx_lo_1 + j + num_ghosts_1_u_qi[qi])*ghostcell_dim_0_u_qi[qi] +
                                            (relative_idx_lo_2 + k + num_ghosts_2_u_qi[qi])*ghostcell_dim_0_u_qi[qi]*
                                                ghostcell_dim_1_u_qi[qi];
                                        
                                        if (use_reciprocal[qi])
                                        {
                                            corr *= (Real(1)/(u_qi[qi][idx_qi]) - u_qi_avg_global[qi][idx_fine]);
                                        }
                                        else
                                        {
                                            corr *= (u_qi[qi][idx_qi] - u_qi_avg_global[qi][idx_fine]);
                                        }
                                    }
                                    
                                    corr_local[idx_fine] += (corr*weight/((Real) n_overlapped));
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                unregisterPatch();
            }
        }
        
        d_mpi.Allreduce(
            corr_local,
            corr_global,
            finest_level_dim_2,
            HAMERS_MPI_REAL,
            MPI_SUM);
        
        std::free(corr_local);
    }
    
    return correlation;
}
