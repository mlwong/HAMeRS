#include "util/MPI_helpers/MPIHelperNumberOfCells.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

/*
 * Compute number of cells.
 */
double
MPIHelperNumberOfCells::getNumberOfCells() const
{
    double num_cells_global;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    if (d_dim == tbox::Dimension(1))
    {
        double num_cells_local = double(0);
        num_cells_global       = double(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Add the number of cells in current patch.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector interior_dims = patch_box.numberCells();
                
                num_cells_local += double(interior_dims[0]);
            }
        }
        
        /*
         * Reduction to get the global number of cells.
         */
        
        d_mpi.Allreduce(
            &num_cells_local,
            &num_cells_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double num_cells_local = double(0);
        num_cells_global       = double(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Add the number of cells in current patch.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector interior_dims = patch_box.numberCells();
                
                num_cells_local += double(interior_dims[0])*double(interior_dims[1]);
            }
        }
        
        /*
         * Reduction to get the global integrals.
         */
        
        d_mpi.Allreduce(
            &num_cells_local,
            &num_cells_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double num_cells_local = double(0);
        num_cells_global       = double(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Add the number of cells in current patch.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector interior_dims = patch_box.numberCells();
                
                num_cells_local += double(interior_dims[0])*double(interior_dims[1])*double(interior_dims[2]);
            }
        }
        
        /*
         * Reduction to get the global number of cells.
         */
        
        d_mpi.Allreduce(
            &num_cells_local,
            &num_cells_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    
    return num_cells_global;
}


/*
 * Compute weighted number of cells.
 */
double
MPIHelperNumberOfCells::getWeightedNumberOfCells() const
{
    double weighted_num_cells_global;
    
    const int num_levels = d_patch_hierarchy->getNumberOfLevels();
    
    if (d_dim == tbox::Dimension(1))
    {
        double weighted_num_cells_local = double(0);
        weighted_num_cells_global       = double(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the refinement ratio from the current level to the coarsest level.
             */
            
            hier::IntVector ratioCurrentLevelToCoarsestLevel =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioCurrentLevelToCoarsestLevel *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            if (li == 0)
            {
                ratioCurrentLevelToCoarsestLevel = hier::IntVector::getOne(d_dim);
            }
            
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Add the number of cells in current patch.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector interior_dims = patch_box.numberCells();
                
                weighted_num_cells_local += double(interior_dims[0])*double(ratioCurrentLevelToCoarsestLevel[0]);
            }
        }
        
        /*
         * Reduction to get the global number of cells.
         */
        
        d_mpi.Allreduce(
            &weighted_num_cells_local,
            &weighted_num_cells_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double weighted_num_cells_local = double(0);
        weighted_num_cells_global       = double(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the refinement ratio from the current level to the coarsest level.
             */
            
            hier::IntVector ratioCurrentLevelToCoarsestLevel =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioCurrentLevelToCoarsestLevel *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            if (li == 0)
            {
                ratioCurrentLevelToCoarsestLevel = hier::IntVector::getOne(d_dim);
            }
            
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Add the number of cells in current patch.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector interior_dims = patch_box.numberCells();
                
                weighted_num_cells_local += double(interior_dims[0])*double(interior_dims[1])*
                    double(ratioCurrentLevelToCoarsestLevel[0]);
            }
        }
        
        /*
         * Reduction to get the global integrals.
         */
        
        d_mpi.Allreduce(
            &weighted_num_cells_local,
            &weighted_num_cells_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double weighted_num_cells_local = double(0);
        weighted_num_cells_global       = double(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the refinement ratio from the current level to the coarsest level.
             */
            
            hier::IntVector ratioCurrentLevelToCoarsestLevel =
                d_patch_hierarchy->getRatioToCoarserLevel(li);
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioCurrentLevelToCoarsestLevel *= d_patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            if (li == 0)
            {
                ratioCurrentLevelToCoarsestLevel = hier::IntVector::getOne(d_dim);
            }
            
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                d_patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Add the number of cells in current patch.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector interior_dims = patch_box.numberCells();
                
                weighted_num_cells_local += double(interior_dims[0])*double(interior_dims[1])*
                    double(ratioCurrentLevelToCoarsestLevel[0]);
            }
        }
        
        /*
         * Reduction to get the global number of cells.
         */
        
        d_mpi.Allreduce(
            &weighted_num_cells_local,
            &weighted_num_cells_global,
            1,
            MPI_DOUBLE,
            MPI_SUM);
    }
    
    return weighted_num_cells_global;
}
