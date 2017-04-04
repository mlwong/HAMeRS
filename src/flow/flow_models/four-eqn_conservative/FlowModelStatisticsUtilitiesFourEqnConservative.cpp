#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "SAMRAI/hier/FlattenedHierarchy.h"

#include <fstream>

/*
 * Output names of statistical quantities to output to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputStatisticalQuantitiesNames(
    const std::string& filename_statistics)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!filename_statistics.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        f_out.open(filename_statistics.c_str(), std::ios::app);
        
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        // Loop over statistical quantities.
        for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
        {
            // Get the key of the current variable.
            std::string statistical_quantity_key = d_statistical_quantities[qi];
            
            if (statistical_quantity_key == "MIXING_WIDTH_X")
            {
                f_out << "\t" << "MIXING_WIDTH_X   ";
            }
            else if (statistical_quantity_key == "MIXING_WIDTH_Y")
            {
                f_out << "\t" << "MIXING_WIDTH_Y   ";
            }
            else if (statistical_quantity_key == "MIXING_WIDTH_Z")
            {
                f_out << "\t" << "MIXING_WIDTH_Z   ";
            }
            else if (statistical_quantity_key == "MIXEDNESS_X")
            {
                f_out << "\t" << "MIXEDNESS_X      ";
            }
            else if (statistical_quantity_key == "MIXEDNESS_Y")
            {
                f_out << "\t" << "MIXEDNESS_Y      ";
            }
            else if (statistical_quantity_key == "MIXEDNESS_Z")
            {
                f_out << "\t" << "MIXEDNESS_Z     ";
            }
            else if (statistical_quantity_key == "ENSTROPHY_INT")
            {
                f_out << "\t" << "ENSTROPHY_INT    ";
            }
            else if (statistical_quantity_key == "TKE_INT_HOMO_X")
            {
                f_out << "\t" << "TKE_INT_HOMO_X   ";
            }
            else if (statistical_quantity_key == "TKE_INT_HOMO_Y")
            {
                f_out << "\t" << "TKE_INT_HOMO_Y   ";
            }
            else if (statistical_quantity_key == "TKE_INT_HOMO_Z")
            {
                f_out << "\t" << "TKE_INT_HOMO_Z   ";
            }
        }
        
        f_out.close();
    }
}


/*
 * Output statisitcal quantities to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputStatisticalQuantities(
    const std::string& filename_statistics,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!filename_statistics.empty());
#endif
    
    // Loop over statistical quantities.
    for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_statistical_quantities[qi];
        
        if (statistical_quantity_key == "MIXING_WIDTH_X")
        {
            outputMixingWidthInXDirection(
                filename_statistics,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXING_WIDTH_Y")
        {
            outputMixingWidthInYDirection(
                filename_statistics,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXING_WIDTH_Z")
        {
            outputMixingWidthInZDirection(
                filename_statistics,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_X")
        {
            outputMixednessInXDirection(
                filename_statistics,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_Y")
        {
            outputMixednessInYDirection(
                filename_statistics,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_Z")
        {
            outputMixednessInZDirection(
                filename_statistics,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "ENSTROPHY_INT")
        {
            outputEnstrophyIntegrated(
                filename_statistics,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_INT_HOMO_X")
        {
            outputTKEIntegratedWithHomogeneityInXDirection(
                filename_statistics,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_INT_HOMO_Y")
        {
            outputTKEIntegratedWithHomogeneityInYDirection(
                filename_statistics,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_INT_HOMO_Z")
        {
            outputTKEIntegratedWithHomogeneityInZDirection(
                filename_statistics,
                patch_hierarchy,
                data_context);
        }
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        f_out.open(filename_statistics.c_str(), std::ios::app);
        
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
        
        f_out.close();
    }
}


/*
 * Output mixing width in x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputMixingWidthInXDirection(
    const std::string& filename_statistics,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!filename_statistics.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(filename_statistics.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<hier::FlattenedHierarchy> flattened_hierarchy(
        new hier::FlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXING_WIDTH_X' can be computed with two species only."
            << std::endl);
    }
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_x = x_hi[0] - x_lo[0];
        
        double* Y_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* Y_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            Y_avg_local[i] = 0.0;
            Y_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and mass fraction in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                double* Y = data_mass_fraction->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const int num_ghosts_0 = num_ghosts[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0;
                        
                        const double value_to_add = Y[idx];
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            Y_avg_local[idx_fine] += value_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Reduce(
            Y_avg_local,
            Y_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the mixing width (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double W = 0.0;
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:W)
#endif
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                W += Y_avg_global[i]*(1.0 - Y_avg_global[i]);
            }
            
            const double dx_finest = L_x/finest_level_dim_0;
            
            W = 4.0*W*dx_finest;
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << W;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_x = x_hi[0] - x_lo[0];
        const double L_y = x_hi[1] - x_lo[1];
        
        double* Y_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* Y_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            Y_avg_local[i] = 0.0;
            Y_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and mass fraction in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                double* Y = data_mass_fraction->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                
                const double weight = dx[1]/L_y;
                
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
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const double value_to_add = Y[idx]*weight;
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                Y_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Reduce(
            Y_avg_local,
            Y_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the mixing width (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double W = 0.0;
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:W)
#endif
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                W += Y_avg_global[i]*(1.0 - Y_avg_global[i]);
            }
            
            const double dx_finest = L_x/finest_level_dim_0;
            
            W = 4.0*W*dx_finest;
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << W;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_x = x_hi[0] - x_lo[0];
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* Y_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* Y_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            Y_avg_local[i] = 0.0;
            Y_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and mass fraction in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                double* Y = data_mass_fraction->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int num_ghosts_2 = num_ghosts[2];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                const int ghostcell_dim_1 = ghostcell_dims[1];
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
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
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                    (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0 +
                                    (relative_idx_lo_2 + k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                const double value_to_add = Y[idx]*weight;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    Y_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Reduce(
            Y_avg_local,
            Y_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the mixing width (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double W = 0.0;
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:W)
#endif
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                W += Y_avg_global[i]*(1.0 - Y_avg_global[i]);
            }
            
            const double dx_finest = L_x/finest_level_dim_0;
            
            W = 4.0*W*dx_finest;
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << W;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
    }
}


/*
 * Output mixing width in y-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputMixingWidthInYDirection(
    const std::string& filename_statistics,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!filename_statistics.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(filename_statistics.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<hier::FlattenedHierarchy> flattened_hierarchy(
        new hier::FlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXING_WIDTH_Y' can be computed with two species only."
            << std::endl);
    }
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'MIXING_WIDTH_Y' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_1 = finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_x = x_hi[0] - x_lo[0];
        const double L_y = x_hi[1] - x_lo[1];
        
        double* Y_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* Y_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            Y_avg_local[j] = 0.0;
            Y_avg_global[j] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_1 = ratioToFinestLevel[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and mass fraction in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                double* Y = data_mass_fraction->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                
                const double weight = dx[0]/L_x;
                
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
                    
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const double value_to_add = Y[idx]*weight;
                            
                            for (int jj = 0; jj < ratioToFinestLevel_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratioToFinestLevel_1 + jj;
                                
                                Y_avg_local[idx_fine] += value_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Reduce(
            Y_avg_local,
            Y_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the mixing width (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double W = 0.0;
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:W)
#endif
            for (int j = 0; j < finest_level_dim_1; j++)
            {
                W += Y_avg_global[j]*(1.0 - Y_avg_global[j]);
            }
            
            const double dy_finest = L_y/finest_level_dim_1;
            
            W = 4.0*W*dy_finest;
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << W;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_x = x_hi[0] - x_lo[0];
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* Y_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* Y_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            Y_avg_local[j] = 0.0;
            Y_avg_global[j] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_1 = ratioToFinestLevel[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and mass fraction in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                double* Y = data_mass_fraction->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int num_ghosts_2 = num_ghosts[2];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                const int ghostcell_dim_1 = ghostcell_dims[1];
                
                const double weight = (dx[0]*dx[2])/(L_x*L_z);
                
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
                    
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                    (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0 +
                                    (relative_idx_lo_2 + k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                const double value_to_add = Y[idx]*weight;
                                
                                for (int jj = 0; jj < ratioToFinestLevel_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratioToFinestLevel_1 + jj;
                                    
                                    Y_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Reduce(
            Y_avg_local,
            Y_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the mixing width (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double W = 0.0;
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:W)
#endif
            for (int j = 0; j < finest_level_dim_1; j++)
            {
                W += Y_avg_global[j]*(1.0 - Y_avg_global[j]);
            }
            
            const double dy_finest = L_y/finest_level_dim_1;
            
            W = 4.0*W*dy_finest;
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << W;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
    }
}


/*
 * Output mixing width in z-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputMixingWidthInZDirection(
    const std::string& filename_statistics,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!filename_statistics.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(filename_statistics.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<hier::FlattenedHierarchy> flattened_hierarchy(
        new hier::FlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXING_WIDTH_Z' can be computed with two species only."
            << std::endl);
    }
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'MIXING_WIDTH_Z' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'MIXING_WIDTH_Z' for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_2 = finest_level_dims[2];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_x = x_hi[0] - x_lo[0];
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* Y_avg_local = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        double* Y_avg_global = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            Y_avg_local[k] = 0.0;
            Y_avg_global[k] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_2 = ratioToFinestLevel[2];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and mass fraction in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                double* Y = data_mass_fraction->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int num_ghosts_2 = num_ghosts[2];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                const int ghostcell_dim_1 = ghostcell_dims[1];
                
                const double weight = (dx[0]*dx[1])/(L_x*L_y);
                
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
                    
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                    (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0 +
                                    (relative_idx_lo_2 + k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                const double value_to_add = Y[idx]*weight;
                                
                                for (int kk = 0; kk < ratioToFinestLevel_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratioToFinestLevel_2 + kk;
                                    
                                    Y_avg_local[idx_fine] += value_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global average.
         */
        
        mpi.Reduce(
            Y_avg_local,
            Y_avg_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the mixing width (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double W = 0.0;
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:W)
#endif
            for (int k = 0; k < finest_level_dim_2; k++)
            {
                W += Y_avg_global[k]*(1.0 - Y_avg_global[k]);
            }
            
            const double dz_finest = L_z/finest_level_dim_2;
            
            W = 4.0*W*dz_finest;
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << W;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
    }
}


/*
 * Output mixedness in x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputMixednessInXDirection(
    const std::string& filename_statistics,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!filename_statistics.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(filename_statistics.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<hier::FlattenedHierarchy> flattened_hierarchy(
        new hier::FlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXEDNESS_X' can be computed with two species only."
            << std::endl);
    }
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        double* Y_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* Y_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* Y_product_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* Y_product_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            Y_avg_local[i] = 0.0;
            Y_avg_global[i] = 0.0;
            Y_product_avg_local[i] = 0.0;
            Y_product_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and mass fraction in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                double* Y_0 = data_mass_fraction->getPointer(0);
                double* Y_1 = data_mass_fraction->getPointer(1);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const int num_ghosts_0 = num_ghosts[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int idx_lo_0 = index_lo[0];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0;
                        
                        const double value_to_add = Y_0[idx];
                        const double product_to_add = Y_0[idx]*Y_1[idx];
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            Y_avg_local[idx_fine] += value_to_add;
                            Y_product_avg_local[idx_fine] += product_to_add;
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global averages.
         */
        
        mpi.Reduce(
            Y_avg_local,
            Y_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            Y_product_avg_local,
            Y_product_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the mixing width (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double num = 0.0;
            double den = 0.0;
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:num)
#endif
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                num += Y_product_avg_global[i];
            }
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:den)
#endif
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                den += Y_avg_global[i]*(1.0 - Y_avg_global[i]);
            }
            
            const double Theta = num/den;
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Theta;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
        std::free(Y_product_avg_local);
        std::free(Y_product_avg_global);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* Y_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* Y_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* Y_product_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* Y_product_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            Y_avg_local[i] = 0.0;
            Y_avg_global[i] = 0.0;
            Y_product_avg_local[i] = 0.0;
            Y_product_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and mass fraction in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                double* Y_0 = data_mass_fraction->getPointer(0);
                double* Y_1 = data_mass_fraction->getPointer(1);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                
                const double weight = dx[1]/L_y;
                
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
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const double value_to_add = Y_0[idx]*weight;
                            const double product_to_add = Y_0[idx]*Y_1[idx]*weight;
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                Y_avg_local[idx_fine] += value_to_add;
                                Y_product_avg_local[idx_fine] += product_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global averages.
         */
        
        mpi.Reduce(
            Y_avg_local,
            Y_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            Y_product_avg_local,
            Y_product_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the mixing width (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double num = 0.0;
            double den = 0.0;
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:num)
#endif
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                num += Y_product_avg_global[i];
            }
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:den)
#endif
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                den += Y_avg_global[i]*(1.0 - Y_avg_global[i]);
            }
            
            const double Theta = num/den;
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Theta;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
        std::free(Y_product_avg_local);
        std::free(Y_product_avg_global);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* Y_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* Y_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* Y_product_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* Y_product_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            Y_avg_local[i] = 0.0;
            Y_avg_global[i] = 0.0;
            Y_product_avg_local[i] = 0.0;
            Y_product_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and mass fraction in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                double* Y_0 = data_mass_fraction->getPointer(0);
                double* Y_1 = data_mass_fraction->getPointer(1);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int num_ghosts_2 = num_ghosts[2];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                const int ghostcell_dim_1 = ghostcell_dims[1];
                
                const double weight = (dx[1]*dx[2])/(L_y*L_z);
                
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
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                    (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0 +
                                    (relative_idx_lo_2 + k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                const double value_to_add = Y_0[idx]*weight;
                                const double product_to_add = Y_0[idx]*Y_1[idx]*weight;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    Y_avg_local[idx_fine] += value_to_add;
                                    Y_product_avg_local[idx_fine] += product_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global averages.
         */
        
        mpi.Reduce(
            Y_avg_local,
            Y_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            Y_product_avg_local,
            Y_product_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the mixing width (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double num = 0.0;
            double den = 0.0;
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:num)
#endif
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                num += Y_product_avg_global[i];
            }
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:den)
#endif
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                den += Y_avg_global[i]*(1.0 - Y_avg_global[i]);
            }
            
            const double Theta = num/den;
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Theta;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
        std::free(Y_product_avg_local);
        std::free(Y_product_avg_global);
    }
}


/*
 * Output mixedness in y-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputMixednessInYDirection(
    const std::string& filename_statistics,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!filename_statistics.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(filename_statistics.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<hier::FlattenedHierarchy> flattened_hierarchy(
        new hier::FlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXEDNESS_Y' can be computed with two species only."
            << std::endl);
    }
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'MIXEDNESS_Y' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_1 = finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_x = x_hi[0] - x_lo[0];
        
        double* Y_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* Y_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* Y_product_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* Y_product_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            Y_avg_local[j] = 0.0;
            Y_avg_global[j] = 0.0;
            Y_product_avg_local[j] = 0.0;
            Y_product_avg_global[j] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_1 = ratioToFinestLevel[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and mass fraction in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                double* Y_0 = data_mass_fraction->getPointer(0);
                double* Y_1 = data_mass_fraction->getPointer(1);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                
                const double weight = dx[0]/L_x;
                
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
                    
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const double value_to_add = Y_0[idx]*weight;
                            const double product_to_add = Y_0[idx]*Y_1[idx]*weight;
                            
                            for (int jj = 0; jj < ratioToFinestLevel_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratioToFinestLevel_1 + jj;
                                
                                Y_avg_local[idx_fine] += value_to_add;
                                Y_product_avg_local[idx_fine] += product_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global averages.
         */
        
        mpi.Reduce(
            Y_avg_local,
            Y_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            Y_product_avg_local,
            Y_product_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the mixing width (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double num = 0.0;
            double den = 0.0;
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:num)
#endif
            for (int j = 0; j < finest_level_dim_1; j++)
            {
                num += Y_product_avg_global[j];
            }
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:den)
#endif
            for (int j = 0; j < finest_level_dim_1; j++)
            {
                den += Y_avg_global[j]*(1.0 - Y_avg_global[j]);
            }
            
            const double Theta = num/den;
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Theta;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
        std::free(Y_product_avg_local);
        std::free(Y_product_avg_global);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_x = x_hi[0] - x_lo[0];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* Y_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* Y_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* Y_product_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* Y_product_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            Y_avg_local[j] = 0.0;
            Y_avg_global[j] = 0.0;
            Y_product_avg_local[j] = 0.0;
            Y_product_avg_global[j] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_1 = ratioToFinestLevel[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and mass fraction in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                double* Y_0 = data_mass_fraction->getPointer(0);
                double* Y_1 = data_mass_fraction->getPointer(1);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int num_ghosts_2 = num_ghosts[2];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                const int ghostcell_dim_1 = ghostcell_dims[1];
                
                const double weight = (dx[0]*dx[2])/(L_x*L_z);
                
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
                    
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                    (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0 +
                                    (relative_idx_lo_2 + k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                const double value_to_add = Y_0[idx]*weight;
                                const double product_to_add = Y_0[idx]*Y_1[idx]*weight;
                                
                                for (int jj = 0; jj < ratioToFinestLevel_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratioToFinestLevel_1 + jj;
                                    
                                    Y_avg_local[idx_fine] += value_to_add;
                                    Y_product_avg_local[idx_fine] += product_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global averages.
         */
        
        mpi.Reduce(
            Y_avg_local,
            Y_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            Y_product_avg_local,
            Y_product_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the mixing width (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double num = 0.0;
            double den = 0.0;
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:num)
#endif
            for (int j = 0; j < finest_level_dim_1; j++)
            {
                num += Y_product_avg_global[j];
            }
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:den)
#endif
            for (int j = 0; j < finest_level_dim_1; j++)
            {
                den += Y_avg_global[j]*(1.0 - Y_avg_global[j]);
            }
            
            const double Theta = num/den;
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Theta;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
        std::free(Y_product_avg_local);
        std::free(Y_product_avg_global);
    }
}


/*
 * Output mixedness in z-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputMixednessInZDirection(
    const std::string& filename_statistics,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!filename_statistics.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(filename_statistics.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<hier::FlattenedHierarchy> flattened_hierarchy(
        new hier::FlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXEDNESS_Z' can be computed with two species only."
            << std::endl);
    }
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'MIXEDNESS_Z' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'MIXEDNESS_Z' for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_2 = finest_level_dims[2];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_x = x_hi[0] - x_lo[0];
        const double L_y = x_hi[1] - x_lo[1];
        
        double* Y_avg_local = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        double* Y_avg_global = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        double* Y_product_avg_local = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        double* Y_product_avg_global = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            Y_avg_local[k] = 0.0;
            Y_avg_global[k] = 0.0;
            Y_product_avg_local[k] = 0.0;
            Y_product_avg_global[k] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_2 = ratioToFinestLevel[2];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and mass fraction in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointer to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                double* Y_0 = data_mass_fraction->getPointer(0);
                double* Y_1 = data_mass_fraction->getPointer(1);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int num_ghosts_2 = num_ghosts[2];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                const int ghostcell_dim_1 = ghostcell_dims[1];
                
                const double weight = (dx[0]*dx[1])/(L_x*L_y);
                
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
                    
                    const int idx_lo_2 = index_lo[2];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                    (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0 +
                                    (relative_idx_lo_2 + k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                const double value_to_add = Y_0[idx]*weight;
                                const double product_to_add = Y_0[idx]*Y_1[idx]*weight;
                                
                                for (int kk = 0; kk < ratioToFinestLevel_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratioToFinestLevel_2 + kk;
                                    
                                    Y_avg_local[idx_fine] += value_to_add;
                                    Y_product_avg_local[idx_fine] += product_to_add;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global averages.
         */
        
        mpi.Reduce(
            Y_avg_local,
            Y_avg_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            Y_product_avg_local,
            Y_product_avg_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the mixing width (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            double num = 0.0;
            double den = 0.0;
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:num)
#endif
            for (int k = 0; k < finest_level_dim_2; k++)
            {
                num += Y_product_avg_global[k];
            }
            
#ifdef HAMERS_ENABLE_SIMD
            #pragma omp simd reduction(+:den)
#endif
            for (int k = 0; k < finest_level_dim_2; k++)
            {
                den += Y_avg_global[k]*(1.0 - Y_avg_global[k]);
            }
            
            const double Theta = num/den;
            
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Theta;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
        std::free(Y_product_avg_local);
        std::free(Y_product_avg_global);
    }
}


/*
 * Output enstrophy integrated to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputEnstrophyIntegrated(
    const std::string& filename_statistics,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!filename_statistics.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(filename_statistics.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<hier::FlattenedHierarchy> flattened_hierarchy(
        new hier::FlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'ENSTROPHY_INT' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double Omega_integrated_local = 0.0;
        double Omega_integrated_global = 0.0;
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch geometry.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch, density and enstrophy in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("ENSTROPHY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to density and enstrophy data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_enstrophy =
                    d_flow_model_tmp->getGlobalCellData("ENSTROPHY");
                
                double* rho = data_density->getPointer(0);
                double* Omega = data_enstrophy->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                
                double Omega_to_add = 0.0;
                
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
                    
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd reduction(+:Omega_to_add)
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0;
                            
                            Omega_to_add += rho[idx]*Omega[idx];
                        }
                    }
                }
                
                Omega_to_add = Omega_to_add*dx[0]*dx[1];
                Omega_integrated_local += Omega_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global integral.
         */
        
        mpi.Reduce(
            &Omega_integrated_local,
            &Omega_integrated_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the enstrophy integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Omega_integrated_global;
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double Omega_integrated_local = 0.0;
        double Omega_integrated_global = 0.0;
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch geometry.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch, density and enstrophy in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("ENSTROPHY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to density and enstrophy data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_enstrophy =
                    d_flow_model_tmp->getGlobalCellData("ENSTROPHY");
                
                double* rho = data_density->getPointer(0);
                double* Omega = data_enstrophy->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int num_ghosts_2 = num_ghosts[2];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                const int ghostcell_dim_1 = ghostcell_dims[1];
                
                double Omega_to_add = 0.0;
                
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
                    
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    const int relative_idx_lo_2 = relative_index_lo[2];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd reduction(+:Omega_to_add)
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                    (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0 +
                                    (relative_idx_lo_2 + k + num_ghosts_2)*ghostcell_dim_0*
                                        ghostcell_dim_1;
                                
                                Omega_to_add += rho[idx]*Omega[idx];
                            }
                        }
                    }
                }
                
                Omega_to_add = Omega_to_add*dx[0]*dx[1]*dx[2];
                Omega_integrated_local += Omega_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global integral.
         */
        
        mpi.Reduce(
            &Omega_integrated_local,
            &Omega_integrated_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the enstrophy integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Omega_integrated_global;
        }
    }
}


/*
 * Output TKE integrated with assumed homogeneity in x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputTKEIntegratedWithHomogeneityInXDirection(
    const std::string& filename_statistics,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!filename_statistics.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(filename_statistics.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<hier::FlattenedHierarchy> flattened_hierarchy(
        new hier::FlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'TKE_INT_HOMO_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_1 = finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_x = x_hi[0] - x_lo[0];
        
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* rho_u_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* rho_u_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* rho_v_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* rho_v_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            rho_avg_local[j] = 0.0;
            rho_avg_global[j] = 0.0;
            rho_u_avg_local[j] = 0.0;
            rho_u_avg_global[j] = 0.0;
            rho_v_avg_local[j] = 0.0;
            rho_v_avg_global[j] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_1 = ratioToFinestLevel[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch, density and velocity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                
                const double weight = dx[0]/L_x;
                
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
                    
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const double rho_to_add = rho[idx]*weight;
                            const double rho_u_to_add = rho[idx]*u[idx]*weight;
                            const double rho_v_to_add = rho[idx]*v[idx]*weight;
                            
                            for (int jj = 0; jj < ratioToFinestLevel_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratioToFinestLevel_1 + jj;
                                
                                rho_avg_local[idx_fine] += rho_to_add;
                                rho_u_avg_local[idx_fine] += rho_u_to_add;
                                rho_v_avg_local[idx_fine] += rho_v_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_u_avg_local,
            rho_u_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_v_avg_local,
            rho_v_avg_global,
            finest_level_dim_1,
            MPI_DOUBLE,
            MPI_SUM);
        
        double TKE_integrated_local = 0.0;
        double TKE_integrated_global = 0.0;
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_1 = ratioToFinestLevel[1];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch, density and velocity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                
                double TKE_to_add = 0.0;
                
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
                    
                    const int idx_lo_1 = index_lo[1];
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd reduction(+:TKE_to_add)
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0;
                            
                            for (int jj = 0; jj < ratioToFinestLevel_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratioToFinestLevel_1 + jj;
                                
                                double u_prime = u[idx] - rho_u_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                double v_prime = v[idx] - rho_v_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                TKE_to_add += 0.5*rho[idx]*(u_prime*u_prime + v_prime*v_prime);
                            }
                        }
                    }
                }
                
                TKE_to_add = TKE_to_add*dx[0]*dx[1]/ratioToFinestLevel_1;
                TKE_integrated_local += TKE_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        mpi.Reduce(
            &TKE_integrated_local,
            &TKE_integrated_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        if (mpi.getRank() == 0)
        {
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << TKE_integrated_global;
        }
        
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_u_avg_local);
        std::free(rho_u_avg_global);
        std::free(rho_v_avg_local);
        std::free(rho_v_avg_global);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'TKE_INT_HOMO_X' is not implemented for three-dimensional problem."
            << std::endl);
    }
}


/*
 * Output TKE integrated with assumed homogeneity in y-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputTKEIntegratedWithHomogeneityInYDirection(
    const std::string& filename_statistics,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!filename_statistics.empty());
#endif
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> d_flow_model_tmp = d_flow_model.lock();
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    
    std::ofstream f_out;
    
    if (mpi.getRank() == 0)
    {
        f_out.open(filename_statistics.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    /*
     * Get the refinement ratio from the finest level to the coarest level.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    hier::IntVector ratioFinestLevelToCoarestLevel =
        patch_hierarchy->getRatioToCoarserLevel(num_levels - 1);
    for (int li = num_levels - 2; li > 0 ; li--)
    {
        ratioFinestLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(li);
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    boost::shared_ptr<hier::FlattenedHierarchy> flattened_hierarchy(
        new hier::FlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the number of cells of physical domain refined to the finest level.
     */
    
    const hier::BoxContainer& physical_domain = d_grid_geometry->getPhysicalDomain();
    const hier::Box& physical_domain_box = physical_domain.front();
    const hier::IntVector& physical_domain_dims = physical_domain_box.numberCells();
    const hier::IntVector finest_level_dims = physical_domain_dims*ratioFinestLevelToCoarestLevel;
    
    /*
     * Get the indices of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'TKE_INT_HOMO_Y' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
            rho_u_avg_local[i] = 0.0;
            rho_u_avg_global[i] = 0.0;
            rho_v_avg_local[i] = 0.0;
            rho_v_avg_global[i] = 0.0;
        }
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch, density and velocity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                
                const double weight = dx[1]/L_y;
                
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
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0;
                            
                            const double rho_to_add = rho[idx]*weight;
                            const double rho_u_to_add = rho[idx]*u[idx]*weight;
                            const double rho_v_to_add = rho[idx]*v[idx]*weight;
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                rho_avg_local[idx_fine] += rho_to_add;
                                rho_u_avg_local[idx_fine] += rho_u_to_add;
                                rho_v_avg_local[idx_fine] += rho_v_to_add;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_u_avg_local,
            rho_u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_v_avg_local,
            rho_v_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        double TKE_integrated_local = 0.0;
        double TKE_integrated_global = 0.0;
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            /*
             * Get the refinement ratio from current level to the finest level.
             */
            
            hier::IntVector ratioToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            hier::IntVector ratioToFinestLevel = ratioFinestLevelToCoarestLevel/ratioToCoarestLevel;
            
            const int ratioToFinestLevel_0 = ratioToFinestLevel[0];
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
                // Get the dimensions of box that covers the interior of patch.
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
                    BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch, density and velocity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    d_flow_model_tmp->getGlobalCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    d_flow_model_tmp->getGlobalCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                double* v = data_velocity->getPointer(1);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                hier::Box ghost_box = patch_box;
                ghost_box.grow(num_ghosts);
                const hier::IntVector ghostcell_dims = ghost_box.numberCells();
                
                const int num_ghosts_0 = num_ghosts[0];
                const int num_ghosts_1 = num_ghosts[1];
                const int ghostcell_dim_0 = ghostcell_dims[0];
                
                double TKE_to_add = 0.0;
                
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
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    const int relative_idx_lo_1 = relative_index_lo[1];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd reduction(+:TKE_to_add)
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0) +
                                (relative_idx_lo_1 + j + num_ghosts_1)*ghostcell_dim_0;
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                double u_prime = u[idx] - rho_u_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                double v_prime = v[idx] - rho_v_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                TKE_to_add += 0.5*rho[idx]*(u_prime*u_prime + v_prime*v_prime);
                            }
                        }
                    }
                }
                
                TKE_to_add = TKE_to_add*dx[0]*dx[1]/ratioToFinestLevel_0;
                TKE_integrated_local += TKE_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                d_flow_model_tmp->unregisterPatch();
            }
        }
        
        mpi.Reduce(
            &TKE_integrated_local,
            &TKE_integrated_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        if (mpi.getRank() == 0)
        {
            f_out << std::fixed << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << TKE_integrated_global;
        }
        
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_u_avg_local);
        std::free(rho_u_avg_global);
        std::free(rho_v_avg_local);
        std::free(rho_v_avg_global);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'TKE_INT_HOMO_Y' is not implemented for three-dimensional problem."
            << std::endl);
    }
}


/*
 * Output TKE integrated with assumed homogeneity in z-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputTKEIntegratedWithHomogeneityInZDirection(
    const std::string& filename_statistics,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'TKE_INT_HOMO_Z' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'TKE_INT_HOMO_Z' for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'TKE_INT_HOMO_Z' is not implemented for three-dimensional problem."
            << std::endl);
    }
}
