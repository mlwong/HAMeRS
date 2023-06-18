#include "flow/refinement_taggers/ImmersedBoundaryTagger.hpp"

#include "util/immersed_boundaries/ImmersedBoundaries.hpp"

ImmersedBoundaryTagger::ImmersedBoundaryTagger(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const HAMERS_SHARED_PTR<FlowModel>& flow_model,
    const hier::IntVector& num_cells_buffer):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_flow_model(flow_model),
        d_num_cells_buffer(num_cells_buffer),
        d_num_immersed_boundary_tagger_ghosts(d_num_cells_buffer)
{
}


/*
 * Print all characteristics of the immersed boundary tagger class.
 */
void
ImmersedBoundaryTagger::printClassData(std::ostream& os) const
{
    os << "\nPrint GradientTagger object..."
       << std::endl;
    
    os << std::endl;
    
    os << "d_num_cells_buffer = " << d_num_cells_buffer;
}


// /*
//  * Put the characteristics of the immersed boundary tagger class into the restart database.
//  */
// void
// ImmersedBoundaryTagger::putToRestart(
//     const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
// {
//     NULL_USE(restart_db);
// }


/*
 * Tag cells on a patch for refinement using immersed boundary masks.
 */
void
ImmersedBoundaryTagger::tagCellsOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(tags->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    d_flow_model->registerPatchWithDataContext(patch, data_context);
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    /*
     * Get the cell data of the immersed boundary mask in the registered patch.
     */
    
    d_flow_model->setupImmersedBoundaryMethod();
    
    HAMERS_SHARED_PTR<FlowModelImmersedBoundaryMethod> flow_model_immersed_boundary_method =
        d_flow_model->getFlowModelImmersedBoundaryMethod();
    
    HAMERS_SHARED_PTR<pdat::CellData<int> > IB_mask_cell_data = flow_model_immersed_boundary_method->
            getCellDataOfImmersedBoundaryMask(data_context);
    
    // Get the number of ghost cells and dimensions of box that covers interior of patch plus
    // ghost cells.
    hier::IntVector num_ghosts_IB_mask     = IB_mask_cell_data->getGhostCellWidth();
    hier::IntVector ghostcell_dims_IB_mask = IB_mask_cell_data->getGhostBox().numberCells();
    
    // Get the pointers to the tag and the immersed boundary mask data.
    int* tag_ptr     = tags->getPointer(0);
    int* IB_mask_ptr = IB_mask_cell_data->getPointer(0);
    
    if (d_dim == tbox::Dimension(1))
    {
        const int interior_dim_0 = interior_dims[0];
        
        const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
        
        HAMERS_PRAGMA_SIMD
        for (int i = 0; i < interior_dim_0; i++)
        {
            // Compute indices.
            const int idx = i + num_ghosts_0_IB_mask;
            const int idx_nghost = i;
            
            if (IB_mask_ptr[idx] == int(IB_MASK::IB_GHOST))
            {
                tag_ptr[idx_nghost] |= 1;
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
        const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
        const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
        
        for (int j = 0; j < interior_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = 0; i < interior_dim_0; i++)
            {
                // Compute indices.
                const int idx = (i + num_ghosts_0_IB_mask) +
                    (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask;
                
                const int idx_nghost = i + j*interior_dim_0;
                
                if (IB_mask_ptr[idx] == int(IB_MASK::IB_GHOST))
                {
                    tag_ptr[idx_nghost] |= 1;
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_ghosts_0_IB_mask = num_ghosts_IB_mask[0];
        const int num_ghosts_1_IB_mask = num_ghosts_IB_mask[1];
        const int num_ghosts_2_IB_mask = num_ghosts_IB_mask[2];
        const int ghostcell_dim_0_IB_mask = ghostcell_dims_IB_mask[0];
        const int ghostcell_dim_1_IB_mask = ghostcell_dims_IB_mask[1];
        
        for (int k = 0; k < interior_dim_2; k++)
        {
            for (int j = 0; j < interior_dim_1; j++)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute indices.
                    const int idx = (i + num_ghosts_0_IB_mask) +
                        (j + num_ghosts_1_IB_mask)*ghostcell_dim_0_IB_mask +
                        (k + num_ghosts_2_IB_mask)*ghostcell_dim_0_IB_mask*
                            ghostcell_dim_1_IB_mask;
                    
                    const int idx_nghost = i + j*interior_dim_0 + k*interior_dim_0*interior_dim_1;
                    
                    if (IB_mask_ptr[idx] == int(IB_MASK::IB_GHOST))
                    {
                        tag_ptr[idx_nghost] |= 1;
                    }
                }
            }
        }
    }
    
    d_flow_model->unregisterPatch();
}
