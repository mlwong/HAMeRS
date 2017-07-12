#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"

#include <fstream>

/*
 * Output names of statistical quantities to output to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputStatisticalQuantitiesNames(
    const std::string& stat_dump_filename)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        
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
                f_out << "\t" << "MIXING_WIDTH_X       ";
            }
            else if (statistical_quantity_key == "MIXING_WIDTH_Y")
            {
                f_out << "\t" << "MIXING_WIDTH_Y       ";
            }
            else if (statistical_quantity_key == "MIXING_WIDTH_Z")
            {
                f_out << "\t" << "MIXING_WIDTH_Z       ";
            }
            else if (statistical_quantity_key == "MIXEDNESS_X")
            {
                f_out << "\t" << "MIXEDNESS_X          ";
            }
            else if (statistical_quantity_key == "MIXEDNESS_Y")
            {
                f_out << "\t" << "MIXEDNESS_Y          ";
            }
            else if (statistical_quantity_key == "MIXEDNESS_Z")
            {
                f_out << "\t" << "MIXEDNESS_Z          ";
            }
            else if (statistical_quantity_key == "TKE_INT_HOMO_X")
            {
                f_out << "\t" << "TKE_INT_HOMO_X       ";
            }
            else if (statistical_quantity_key == "TKE_INT_HOMO_Y")
            {
                f_out << "\t" << "TKE_INT_HOMO_Y       ";
            }
            else if (statistical_quantity_key == "TKE_INT_HOMO_Z")
            {
                f_out << "\t" << "TKE_INT_HOMO_Z       ";
            }
            else if (statistical_quantity_key == "TKE_INT_HOMO_XY")
            {
                f_out << "\t" << "TKE_INT_HOMO_XY      ";
            }
            else if (statistical_quantity_key == "TKE_INT_HOMO_YZ")
            {
                f_out << "\t" << "TKE_INT_HOMO_YZ      ";
            }
            else if (statistical_quantity_key == "TKE_INT_HOMO_XZ")
            {
                f_out << "\t" << "TKE_INT_HOMO_XZ      ";
            }
            else if (statistical_quantity_key == "ENSTROPHY_INT")
            {
                f_out << "\t" << "ENSTROPHY_INT        ";
            }
            else if (statistical_quantity_key == "SCAL_DISS_RAT_INT")
            {
                f_out << "\t" << "SCAL_DISS_RAT_INT    ";
            }
            else if (statistical_quantity_key == "NUM_INTEF_THICK")
            {
                f_out << "\t" << "NUM_INTEF_THICK      ";
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
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    // Loop over statistical quantities.
    for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_statistical_quantities[qi];
        
        if (statistical_quantity_key == "MIXING_WIDTH_X")
        {
            outputMixingWidthInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXING_WIDTH_Y")
        {
            outputMixingWidthInYDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXING_WIDTH_Z")
        {
            outputMixingWidthInZDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_X")
        {
            outputMixednessInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_Y")
        {
            outputMixednessInYDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_Z")
        {
            outputMixednessInZDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_INT_HOMO_X")
        {
            outputTKEIntegratedWithHomogeneityInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_INT_HOMO_Y")
        {
            outputTKEIntegratedWithHomogeneityInYDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_INT_HOMO_Z")
        {
            outputTKEIntegratedWithHomogeneityInZDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_INT_HOMO_XY")
        {
            outputTKEIntegratedWithHomogeneityInXYPlane(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_INT_HOMO_YZ")
        {
            outputTKEIntegratedWithHomogeneityInYZPlane(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "TKE_INT_HOMO_XZ")
        {
            outputTKEIntegratedWithHomogeneityInXZPlane(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "ENSTROPHY_INT")
        {
            outputEnstrophyIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "SCAL_DISS_RAT_INT")
        {
            outputScalarDissipationRateIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "NUM_INTEF_THICK")
        {
            outputNumericalInterfaceThickness(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Unknown statistical quantity key = '"
                << statistical_quantity_key
                << " found."
                << std::endl);
        }
    }
    
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    if (mpi.getRank() == 0)
    {
        std::ofstream f_out;
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        
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
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXING_WIDTH_X' can be computed with two species only."
            << std::endl);
    }
    
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                
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
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
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
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_mass_fraction;
                        
                        const double value_to_add = fmax(fmin(Y[idx], 1.0), 0.0)/((double) n_overlapped);
                        
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
            
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                W += Y_avg_global[i]*(1.0 - Y_avg_global[i]);
            }
            
            const double dx_finest = L_x/finest_level_dim_0;
            
            W = 4.0*W*dx_finest;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                
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
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                            
                            const double value_to_add = fmax(fmin(Y[idx], 1.0), 0.0)*weight/((double) n_overlapped);
                            
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
            
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                W += Y_avg_global[i]*(1.0 - Y_avg_global[i]);
            }
            
            const double dx_finest = L_x/finest_level_dim_0;
            
            W = 4.0*W*dx_finest;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int num_ghosts_2_mass_fraction = num_ghosts_mass_fraction[2];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                const int ghostcell_dim_1_mass_fraction = ghostcell_dims_mass_fraction[1];
                
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
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                        ghostcell_dim_1_mass_fraction;
                                
                                const double value_to_add = fmax(fmin(Y[idx], 1.0), 0.0)*weight/((double) n_overlapped);
                                
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
            
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                W += Y_avg_global[i]*(1.0 - Y_avg_global[i]);
            }
            
            const double dx_finest = L_x/finest_level_dim_0;
            
            W = 4.0*W*dx_finest;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXING_WIDTH_Y' can be computed with two species only."
            << std::endl);
    }
    
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                
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
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                            
                            const double value_to_add = fmax(fmin(Y[idx], 1.0), 0.0)*weight/((double) n_overlapped);
                            
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
            
            for (int j = 0; j < finest_level_dim_1; j++)
            {
                W += Y_avg_global[j]*(1.0 - Y_avg_global[j]);
            }
            
            const double dy_finest = L_y/finest_level_dim_1;
            
            W = 4.0*W*dy_finest;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int num_ghosts_2_mass_fraction = num_ghosts_mass_fraction[2];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                const int ghostcell_dim_1_mass_fraction = ghostcell_dims_mass_fraction[1];
                
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
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                        ghostcell_dim_1_mass_fraction;
                                
                                const double value_to_add = fmax(fmin(Y[idx], 1.0), 0.0)*weight/((double) n_overlapped);
                                
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
            
            for (int j = 0; j < finest_level_dim_1; j++)
            {
                W += Y_avg_global[j]*(1.0 - Y_avg_global[j]);
            }
            
            const double dy_finest = L_y/finest_level_dim_1;
            
            W = 4.0*W*dy_finest;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXING_WIDTH_Z' can be computed with two species only."
            << std::endl);
    }
    
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int num_ghosts_2_mass_fraction = num_ghosts_mass_fraction[2];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                const int ghostcell_dim_1_mass_fraction = ghostcell_dims_mass_fraction[1];
                
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
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                        ghostcell_dim_1_mass_fraction;
                                
                                const double value_to_add = fmax(fmin(Y[idx], 1.0), 0.0)*weight/((double) n_overlapped);
                                
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
            
            for (int k = 0; k < finest_level_dim_2; k++)
            {
                W += Y_avg_global[k]*(1.0 - Y_avg_global[k]);
            }
            
            const double dz_finest = L_z/finest_level_dim_2;
            
            W = 4.0*W*dz_finest;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXEDNESS_X' can be computed with two species only."
            << std::endl);
    }
    
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                
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
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the index of the data point and count how many times the data is repeated.
                         */
                        
                        const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                        
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
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_mass_fraction;
                        
                        const double Y_bounded = fmax(fmin(Y[idx], 1.0), 0.0);
                        
                        const double weight_local = 1.0/((double) n_overlapped);
                        
                        const double value_to_add = Y_bounded*weight_local;
                        const double product_to_add = Y_bounded*(1.0 - Y_bounded)*weight_local;
                        
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
            
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                num += Y_product_avg_global[i];
            }
            
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                den += Y_avg_global[i]*(1.0 - Y_avg_global[i]);
            }
            
            const double Theta = num/den;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                
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
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                            
                            const double Y_bounded = fmax(fmin(Y[idx], 1.0), 0.0);
                            
                            const double weight_local = weight/((double) n_overlapped);
                            
                            const double value_to_add = Y_bounded*weight_local;
                            const double product_to_add = Y_bounded*(1.0 - Y_bounded)*weight_local;
                            
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
            
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                num += Y_product_avg_global[i];
            }
            
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                den += Y_avg_global[i]*(1.0 - Y_avg_global[i]);
            }
            
            const double Theta = num/den;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int num_ghosts_2_mass_fraction = num_ghosts_mass_fraction[2];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                const int ghostcell_dim_1_mass_fraction = ghostcell_dims_mass_fraction[1];
                
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
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                        ghostcell_dim_1_mass_fraction;
                                
                                const double Y_bounded = fmax(fmin(Y[idx], 1.0), 0.0);
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double value_to_add = Y_bounded*weight_local;
                                const double product_to_add = Y_bounded*(1.0 - Y_bounded)*weight_local;
                                
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
            
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                num += Y_product_avg_global[i];
            }
            
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                den += Y_avg_global[i]*(1.0 - Y_avg_global[i]);
            }
            
            const double Theta = num/den;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXEDNESS_Y' can be computed with two species only."
            << std::endl);
    }
    
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                
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
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                            
                            const double Y_bounded = fmax(fmin(Y[idx], 1.0), 0.0);
                            
                            const double weight_local = weight/((double) n_overlapped);
                        
                            const double value_to_add = Y_bounded*weight_local;
                            const double product_to_add = Y_bounded*(1.0 - Y_bounded)*weight_local;
                            
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
            
            for (int j = 0; j < finest_level_dim_1; j++)
            {
                num += Y_product_avg_global[j];
            }
            
            for (int j = 0; j < finest_level_dim_1; j++)
            {
                den += Y_avg_global[j]*(1.0 - Y_avg_global[j]);
            }
            
            const double Theta = num/den;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int num_ghosts_2_mass_fraction = num_ghosts_mass_fraction[2];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                const int ghostcell_dim_1_mass_fraction = ghostcell_dims_mass_fraction[1];
                
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
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                        ghostcell_dim_1_mass_fraction;
                                
                                const double Y_bounded = fmax(fmin(Y[idx], 1.0), 0.0);
                                
                                const double weight_local = weight/((double) n_overlapped);
                            
                                const double value_to_add = Y_bounded*weight_local;
                                const double product_to_add = Y_bounded*(1.0 - Y_bounded)*weight_local;
                                
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
            
            for (int j = 0; j < finest_level_dim_1; j++)
            {
                num += Y_product_avg_global[j];
            }
            
            for (int j = 0; j < finest_level_dim_1; j++)
            {
                den += Y_avg_global[j]*(1.0 - Y_avg_global[j]);
            }
            
            const double Theta = num/den;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'MIXEDNESS_Z' can be computed with two species only."
            << std::endl);
    }
    
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int num_ghosts_2_mass_fraction = num_ghosts_mass_fraction[2];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                const int ghostcell_dim_1_mass_fraction = ghostcell_dims_mass_fraction[1];
                
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
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                        ghostcell_dim_1_mass_fraction;
                                
                                const double Y_bounded = fmax(fmin(Y[idx], 1.0), 0.0);
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double value_to_add = Y_bounded*weight_local;
                                const double product_to_add = Y_bounded*(1.0 - Y_bounded)*weight_local;
                                
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
            
            for (int k = 0; k < finest_level_dim_2; k++)
            {
                num += Y_product_avg_global[k];
            }
            
            for (int k = 0; k < finest_level_dim_2; k++)
            {
                den += Y_avg_global[k]*(1.0 - Y_avg_global[k]);
            }
            
            const double Theta = num/den;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Theta;
        }
        
        std::free(Y_avg_local);
        std::free(Y_avg_global);
        std::free(Y_product_avg_local);
        std::free(Y_product_avg_global);
    }
}


/*
 * Output TKE integrated with assumed homogeneity in x-direction to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputTKEIntegratedWithHomogeneityInXDirection(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                
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
                            
                            const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density;
                            
                            const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const double weight_local = weight/((double) n_overlapped);
                            
                            const double rho_to_add = rho[idx_density]*weight_local;
                            const double rho_u_to_add = rho[idx_density]*u[idx_velocity]*weight_local;
                            const double rho_v_to_add = rho[idx_density]*v[idx_velocity]*weight_local;
                            
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
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                
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
                            
                            const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density;
                            
                            const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const double weight = 1.0/((double) n_overlapped);
                            
                            for (int jj = 0; jj < ratioToFinestLevel_1; jj++)
                            {
                                const int idx_fine = (idx_lo_1 + j)*ratioToFinestLevel_1 + jj;
                                
                                double u_prime = u[idx_velocity] -
                                    rho_u_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                
                                double v_prime = v[idx_velocity] -
                                    rho_v_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                
                                TKE_to_add += 0.5*rho[idx_density]*(u_prime*u_prime + v_prime*v_prime)*weight;
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
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                
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
                            
                            const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density;
                            
                            const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const double weight_local = weight/((double) n_overlapped);
                            
                            const double rho_to_add = rho[idx_density]*weight_local;
                            const double rho_u_to_add = rho[idx_density]*u[idx_velocity]*weight_local;
                            const double rho_v_to_add = rho[idx_density]*v[idx_velocity]*weight_local;
                            
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
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                
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
                            
                            const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density;
                            
                            const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity;
                            
                            const double weight = 1.0/((double) n_overlapped);
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                double u_prime = u[idx_velocity] -
                                    rho_u_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                
                                double v_prime = v[idx_velocity] -
                                    rho_v_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                
                                TKE_to_add += 0.5*rho[idx_density]*(u_prime*u_prime + v_prime*v_prime)*weight;
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
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
    const std::string& stat_dump_filename,
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


/*
 * Output TKE integrated with assumed homogeneity in xy-plane to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputTKEIntegratedWithHomogeneityInXYPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
            << "There is no 'TKE_INT_HOMO_XY' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'TKE_INT_HOMO_XY' for two-dimensional problem."
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
        
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        double* rho_u_avg_local = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        double* rho_u_avg_global = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        double* rho_v_avg_local = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        double* rho_v_avg_global = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        double* rho_w_avg_local = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        double* rho_w_avg_global = (double*)std::malloc(finest_level_dim_2*sizeof(double));
        
        for (int k = 0; k < finest_level_dim_2; k++)
        {
            rho_avg_local[k] = 0.0;
            rho_avg_global[k] = 0.0;
            rho_u_avg_local[k] = 0.0;
            rho_u_avg_global[k] = 0.0;
            rho_v_avg_local[k] = 0.0;
            rho_v_avg_global[k] = 0.0;
            rho_w_avg_local[k] = 0.0;
            rho_w_avg_global[k] = 0.0;
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
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
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
                double* w = data_velocity->getPointer(2);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int num_ghosts_2_velocity = num_ghosts_velocity[2];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
                
                const double weight = dx[0]*dx[1]/(L_x*L_y);
                
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
                                 * Compute the linear indices and the data to add.
                                 */
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_u_to_add = rho[idx_density]*u[idx_velocity]*weight_local;
                                const double rho_v_to_add = rho[idx_density]*v[idx_velocity]*weight_local;
                                const double rho_w_to_add = rho[idx_density]*w[idx_velocity]*weight_local;
                                
                                for (int kk = 0; kk < ratioToFinestLevel_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratioToFinestLevel_2 + kk;
                                    
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_u_avg_local[idx_fine] += rho_u_to_add;
                                    rho_v_avg_local[idx_fine] += rho_v_to_add;
                                    rho_w_avg_local[idx_fine] += rho_w_to_add;
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
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_u_avg_local,
            rho_u_avg_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_v_avg_local,
            rho_v_avg_global,
            finest_level_dim_2,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_w_avg_local,
            rho_w_avg_global,
            finest_level_dim_2,
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
            
            const int ratioToFinestLevel_2 = ratioToFinestLevel[2];
            
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
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
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
                double* w = data_velocity->getPointer(2);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int num_ghosts_2_velocity = num_ghosts_velocity[2];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
                
                double TKE_to_add = 0.0;
                
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
                                 * Compute the linear indices and the data to add.
                                 */
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight = 1.0/((double) n_overlapped);
                                
                                for (int kk = 0; kk < ratioToFinestLevel_2; kk++)
                                {
                                    const int idx_fine = (idx_lo_2 + k)*ratioToFinestLevel_2 + kk;
                                    
                                    double u_prime = u[idx_velocity] -
                                        rho_u_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                    
                                    double v_prime = v[idx_velocity] -
                                        rho_v_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                    
                                    double w_prime = w[idx_velocity] -
                                        rho_w_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                    
                                    TKE_to_add += 0.5*rho[idx_density]*
                                        (u_prime*u_prime + v_prime*v_prime + w_prime*w_prime)*weight;
                                }
                            }
                        }
                    }
                }
                
                TKE_to_add = TKE_to_add*dx[0]*dx[1]*dx[2]/ratioToFinestLevel_2;
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
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << TKE_integrated_global;
        }
        
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_u_avg_local);
        std::free(rho_u_avg_global);
        std::free(rho_v_avg_local);
        std::free(rho_v_avg_global);
        std::free(rho_w_avg_local);
        std::free(rho_w_avg_global);
    }
}


/*
 * Output TKE integrated with assumed homogeneity in yz-plane to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputTKEIntegratedWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
            << "There is no 'TKE_INT_HOMO_YZ' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'TKE_INT_HOMO_YZ' for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
            rho_u_avg_local[i] = 0.0;
            rho_u_avg_global[i] = 0.0;
            rho_v_avg_local[i] = 0.0;
            rho_v_avg_global[i] = 0.0;
            rho_w_avg_local[i] = 0.0;
            rho_w_avg_global[i] = 0.0;
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
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
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
                double* w = data_velocity->getPointer(2);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int num_ghosts_2_velocity = num_ghosts_velocity[2];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
                
                const double weight = dx[1]*dx[2]/(L_y*L_z);
                
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
                                 * Compute the linear indices and the data to add.
                                 */
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_u_to_add = rho[idx_density]*u[idx_velocity]*weight_local;
                                const double rho_v_to_add = rho[idx_density]*v[idx_velocity]*weight_local;
                                const double rho_w_to_add = rho[idx_density]*w[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_u_avg_local[idx_fine] += rho_u_to_add;
                                    rho_v_avg_local[idx_fine] += rho_v_to_add;
                                    rho_w_avg_local[idx_fine] += rho_w_to_add;
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
        
        mpi.Allreduce(
            rho_w_avg_local,
            rho_w_avg_global,
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
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
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
                double* w = data_velocity->getPointer(2);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int num_ghosts_2_velocity = num_ghosts_velocity[2];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
                
                double TKE_to_add = 0.0;
                
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
                                 * Compute the linear indices and the data to add.
                                 */
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight = 1.0/((double) n_overlapped);
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    double u_prime = u[idx_velocity] -
                                        rho_u_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                    
                                    double v_prime = v[idx_velocity] -
                                        rho_v_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                    
                                    double w_prime = w[idx_velocity] -
                                        rho_w_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                    
                                    TKE_to_add += 0.5*rho[idx_density]*
                                        (u_prime*u_prime + v_prime*v_prime + w_prime*w_prime)*weight;
                                }
                            }
                        }
                    }
                }
                
                TKE_to_add = TKE_to_add*dx[0]*dx[1]*dx[2]/ratioToFinestLevel_0;
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
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << TKE_integrated_global;
        }
        
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_u_avg_local);
        std::free(rho_u_avg_global);
        std::free(rho_v_avg_local);
        std::free(rho_v_avg_global);
        std::free(rho_w_avg_local);
        std::free(rho_w_avg_global);
    }
}


/*
 * Output TKE integrated with assumed homogeneity in xz-plane to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputTKEIntegratedWithHomogeneityInXZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
            << "There is no 'TKE_INT_HOMO_XZ' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'TKE_INT_HOMO_XZ' for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_1 = finest_level_dims[1];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_x = x_hi[0] - x_lo[0];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* rho_u_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* rho_u_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* rho_v_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* rho_v_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* rho_w_avg_local = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        double* rho_w_avg_global = (double*)std::malloc(finest_level_dim_1*sizeof(double));
        
        for (int j = 0; j < finest_level_dim_1; j++)
        {
            rho_avg_local[j] = 0.0;
            rho_avg_global[j] = 0.0;
            rho_u_avg_local[j] = 0.0;
            rho_u_avg_global[j] = 0.0;
            rho_v_avg_local[j] = 0.0;
            rho_v_avg_global[j] = 0.0;
            rho_w_avg_local[j] = 0.0;
            rho_w_avg_global[j] = 0.0;
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
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
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
                double* w = data_velocity->getPointer(2);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int num_ghosts_2_velocity = num_ghosts_velocity[2];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
                
                const double weight = dx[0]*dx[2]/(L_x*L_z);
                
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
                                 * Compute the linear indices and the data to add.
                                 */
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_u_to_add = rho[idx_density]*u[idx_velocity]*weight_local;
                                const double rho_v_to_add = rho[idx_density]*v[idx_velocity]*weight_local;
                                const double rho_w_to_add = rho[idx_density]*w[idx_velocity]*weight_local;
                                
                                for (int jj = 0; jj < ratioToFinestLevel_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratioToFinestLevel_1 + jj;
                                    
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_u_avg_local[idx_fine] += rho_u_to_add;
                                    rho_v_avg_local[idx_fine] += rho_v_to_add;
                                    rho_w_avg_local[idx_fine] += rho_w_to_add;
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
        
        mpi.Allreduce(
            rho_w_avg_local,
            rho_w_avg_global,
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
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
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
                double* w = data_velocity->getPointer(2);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
                const int num_ghosts_0_velocity = num_ghosts_velocity[0];
                const int num_ghosts_1_velocity = num_ghosts_velocity[1];
                const int num_ghosts_2_velocity = num_ghosts_velocity[2];
                const int ghostcell_dim_0_velocity = ghostcell_dims_velocity[0];
                const int ghostcell_dim_1_velocity = ghostcell_dims_velocity[1];
                
                double TKE_to_add = 0.0;
                
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
                                 * Compute the linear indices and the data to add.
                                 */
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight = 1.0/((double) n_overlapped);
                                
                                for (int jj = 0; jj < ratioToFinestLevel_1; jj++)
                                {
                                    const int idx_fine = (idx_lo_1 + j)*ratioToFinestLevel_1 + jj;
                                    
                                    double u_prime = u[idx_velocity] -
                                        rho_u_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                    
                                    double v_prime = v[idx_velocity] -
                                        rho_v_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                    
                                    double w_prime = w[idx_velocity] -
                                        rho_w_avg_global[idx_fine]/rho_avg_global[idx_fine];
                                    
                                    TKE_to_add += 0.5*rho[idx_density]*
                                        (u_prime*u_prime + v_prime*v_prime + w_prime*w_prime)*weight;
                                }
                            }
                        }
                    }
                }
                
                TKE_to_add = TKE_to_add*dx[0]*dx[1]*dx[2]/ratioToFinestLevel_1;
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
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << TKE_integrated_global;
        }
        
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_u_avg_local);
        std::free(rho_u_avg_global);
        std::free(rho_v_avg_local);
        std::free(rho_v_avg_global);
        std::free(rho_w_avg_local);
        std::free(rho_w_avg_global);
    }
}


/*
 * Output enstrophy integrated to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputEnstrophyIntegrated(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("ENSTROPHY", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_enstrophy = data_enstrophy->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_enstrophy = data_enstrophy->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                
                const int num_ghosts_0_enstrophy = num_ghosts_enstrophy[0];
                const int num_ghosts_1_enstrophy = num_ghosts_enstrophy[1];
                const int ghostcell_dim_0_enstrophy = ghostcell_dims_enstrophy[0];
                
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
                            
                            // Compute linear indices.
                            const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density;
                            
                            const int idx_enstrophy = (relative_idx_lo_0 + i + num_ghosts_0_enstrophy) +
                                (relative_idx_lo_1 + j + num_ghosts_1_enstrophy)*ghostcell_dim_0_enstrophy;
                            
                            Omega_to_add += rho[idx_density]*Omega[idx_enstrophy]/((double) n_overlapped);
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
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
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
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("ENSTROPHY", hier::IntVector::getZero(d_dim)));
                
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
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_enstrophy = data_enstrophy->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_enstrophy = data_enstrophy->getGhostBox().numberCells();
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
                const int num_ghosts_0_enstrophy = num_ghosts_enstrophy[0];
                const int num_ghosts_1_enstrophy = num_ghosts_enstrophy[1];
                const int num_ghosts_2_enstrophy = num_ghosts_enstrophy[2];
                const int ghostcell_dim_0_enstrophy = ghostcell_dims_enstrophy[0];
                const int ghostcell_dim_1_enstrophy = ghostcell_dims_enstrophy[1];
                
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
                                 * Compute the linear indices and the data to add.
                                 */
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_enstrophy = (relative_idx_lo_0 + i + num_ghosts_0_enstrophy) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_enstrophy)*ghostcell_dim_0_enstrophy +
                                    (relative_idx_lo_2 + k + num_ghosts_2_enstrophy)*ghostcell_dim_0_enstrophy*
                                        ghostcell_dim_1_enstrophy;
                                
                                Omega_to_add += rho[idx_density]*Omega[idx_enstrophy]/((double) n_overlapped);
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
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Omega_integrated_global;
        }
    }
}


/*
 * Output scalar dissipation rate of first species integrated to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputScalarDissipationRateIntegrated(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (!d_equation_of_mass_diffusivity_mixing_rules)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Mixing rule of mass diffusivity is not initialized yet."
            << std::endl);
    }
    
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
        if (!f_out.is_open())
        {
            TBOX_ERROR(d_object_name
                << ": "
                << "Failed to open file to output statistics!"
                << std::endl);
        }
    }
    
    /*
     * Get the flattened hierarchy where only the finest existing grid is visible at any given
     * location in the problem space.
     */
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        double Chi_integrated_local = 0.0;
        double Chi_integrated_global = 0.0;
        
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
                 * Register the patch, mass fraction, pressure and temperature in the flow model
                 * and compute the corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to mass fraction, pressure and temperature data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                boost::shared_ptr<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getGlobalCellData("PRESSURE");
                
                boost::shared_ptr<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getGlobalCellData("TEMPERATURE");
                
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                
                double Chi_to_add = 0.0;
                
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
                    
                    if (false) // (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*4)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                            
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
                            
                            // Compute indices of current and neighboring cells.
                            const int idx_mass_fraction = relative_idx_lo_0 + i + num_ghosts_0_mass_fraction;
                            const int idx_x_LLLL = relative_idx_lo_0 + i - 4 + num_ghosts_0_mass_fraction;
                            const int idx_x_LLL = relative_idx_lo_0 + i - 3 + num_ghosts_0_mass_fraction;
                            const int idx_x_LL = relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction;
                            const int idx_x_L = relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_R = relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_RR = relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction;
                            const int idx_x_RRR = relative_idx_lo_0 + i + 3 + num_ghosts_0_mass_fraction;
                            const int idx_x_RRRR = relative_idx_lo_0 + i + 4 + num_ghosts_0_mass_fraction;
                            
                            // Compute linear indices of pressure and temperature.
                            const int idx_pressure = relative_idx_lo_0 + i + num_ghosts_0_pressure;
                            const int idx_temperature = relative_idx_lo_0 + i + num_ghosts_0_temperature;
                            
                            const double dYdx = (-1.0/280.0*Y[0][idx_x_RRRR] + 4.0/105.0*Y[0][idx_x_RRR] -
                                                 1.0/5.0*Y[0][idx_x_RR] + 4.0/5.0*Y[0][idx_x_R] -
                                                 4.0/5.0*Y[0][idx_x_L] + 1.0/5.0*Y[0][idx_x_LL] -
                                                 4.0/105.0*Y[0][idx_x_LLL] + 1.0/280.0*Y[0][idx_x_LLLL])/dx[0];
                            
                            std::vector<double> D;
                            D.resize(d_num_species);
                            
                            std::vector<double*> D_ptr;
                            D_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                D_ptr.push_back(&D[si]);
                            }
                            
                            std::vector<const double*> Y_ptr;
                            Y_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                            }
                            
                            d_equation_of_mass_diffusivity_mixing_rules->
                                getMassDiffusivities(
                                    D_ptr,
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            
                            Chi_to_add += D[0]*dYdx*dYdx/((double) n_overlapped);
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*3)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                            
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
                            
                            // Compute indices of current and neighboring cells.
                            const int idx_mass_fraction = relative_idx_lo_0 + i + num_ghosts_0_mass_fraction;
                            const int idx_x_LLL = relative_idx_lo_0 + i - 3 + num_ghosts_0_mass_fraction;
                            const int idx_x_LL = relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction;
                            const int idx_x_L = relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_R = relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_RR = relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction;
                            const int idx_x_RRR = relative_idx_lo_0 + i + 3 + num_ghosts_0_mass_fraction;
                            
                            // Compute linear indices of pressure and temperature.
                            const int idx_pressure = relative_idx_lo_0 + i + num_ghosts_0_pressure;
                            const int idx_temperature = relative_idx_lo_0 + i + num_ghosts_0_temperature;
                            
                            const double dYdx = (1.0/60.0*Y[0][idx_x_RRR] - 3.0/20.0*Y[0][idx_x_RR] +
                                                 3.0/4.0*Y[0][idx_x_R] - 3.0/4.0*Y[0][idx_x_L] +
                                                 3.0/20.0*Y[0][idx_x_LL] - 1.0/60.0*Y[0][idx_x_LLL])/dx[0];
                            
                            std::vector<double> D;
                            D.resize(d_num_species);
                            
                            std::vector<double*> D_ptr;
                            D_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                D_ptr.push_back(&D[si]);
                            }
                            
                            std::vector<const double*> Y_ptr;
                            Y_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                            }
                            
                            d_equation_of_mass_diffusivity_mixing_rules->
                                getMassDiffusivities(
                                    D_ptr,
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            
                            Chi_to_add += D[0]*dYdx*dYdx/((double) n_overlapped);
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*2)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                            
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
                            
                            // Compute indices of current and neighboring cells.
                            const int idx_mass_fraction = relative_idx_lo_0 + i + num_ghosts_0_mass_fraction;
                            const int idx_x_LL = relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction;
                            const int idx_x_L = relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_R = relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_RR = relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction;
                            
                            // Compute linear indices of pressure and temperature.
                            const int idx_pressure = relative_idx_lo_0 + i + num_ghosts_0_pressure;
                            const int idx_temperature = relative_idx_lo_0 + i + num_ghosts_0_temperature;
                            
                            const double dYdx = (-1.0/12.0*Y[0][idx_x_RR] + 2.0/3.0*Y[0][idx_x_R] -
                                                 2.0/3.0*Y[0][idx_x_L] + 1.0/12.0*Y[0][idx_x_LL])/dx[0];
                            
                            std::vector<double> D;
                            D.resize(d_num_species);
                            
                            std::vector<double*> D_ptr;
                            D_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                D_ptr.push_back(&D[si]);
                            }
                            
                            std::vector<const double*> Y_ptr;
                            Y_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                            }
                            
                            d_equation_of_mass_diffusivity_mixing_rules->
                                getMassDiffusivities(
                                    D_ptr,
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            
                            Chi_to_add += D[0]*dYdx*dYdx/((double) n_overlapped);
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim))
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                            
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
                            
                            // Compute indices of current and neighboring cells.
                            const int idx_mass_fraction = relative_idx_lo_0 + i + num_ghosts_0_mass_fraction;
                            const int idx_x_L = relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_R = relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction;
                            
                            // Compute linear indices of pressure and temperature.
                            const int idx_pressure = relative_idx_lo_0 + i + num_ghosts_0_pressure;
                            const int idx_temperature = relative_idx_lo_0 + i + num_ghosts_0_temperature;
                            
                            const double dYdx = (0.5*Y[0][idx_x_R] - 0.5*Y[0][idx_x_L])/dx[0];
                            
                            std::vector<double> D;
                            D.resize(d_num_species);
                            
                            std::vector<double*> D_ptr;
                            D_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                D_ptr.push_back(&D[si]);
                            }
                            
                            std::vector<const double*> Y_ptr;
                            Y_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                            }
                            
                            d_equation_of_mass_diffusivity_mixing_rules->
                                getMassDiffusivities(
                                    D_ptr,
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            
                            Chi_to_add += D[0]*dYdx*dYdx/((double) n_overlapped);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the index of the data point and count how many times the data is repeated.
                             */
                            
                            const hier::Index idx_pt(tbox::Dimension(1), idx_lo_0 + i);
                            
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
                            
                            // Compute indices of current and neighboring cells.
                            const int idx_mass_fraction = relative_idx_lo_0 + i + num_ghosts_0_mass_fraction;
                            const int idx_x_L = relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_R = relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction;
                            
                            // Compute linear indices of pressure and temperature.
                            const int idx_pressure = relative_idx_lo_0 + i + num_ghosts_0_pressure;
                            const int idx_temperature = relative_idx_lo_0 + i + num_ghosts_0_temperature;
                            
                            double dYdx;
                            
                            if (i == -num_ghosts_0_mass_fraction)
                            {
                                dYdx = (Y[0][idx_x_R] - Y[0][idx_mass_fraction])/dx[0];
                            }
                            else if (i == interior_dim_0 + num_ghosts_0_mass_fraction - 1)
                            {
                                dYdx = (Y[0][idx_mass_fraction] - Y[0][idx_x_L])/dx[0];
                            }
                            else
                            {
                                dYdx = (0.5*Y[0][idx_x_R] - 0.5*Y[0][idx_x_L])/dx[0];
                            }
                            
                            std::vector<double> D;
                            D.resize(d_num_species);
                            
                            std::vector<double*> D_ptr;
                            D_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                D_ptr.push_back(&D[si]);
                            }
                            
                            std::vector<const double*> Y_ptr;
                            Y_ptr.reserve(d_num_species);
                            for (int si = 0; si < d_num_species; si++)
                            {
                                Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                            }
                            
                            d_equation_of_mass_diffusivity_mixing_rules->
                                getMassDiffusivities(
                                    D_ptr,
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            
                            Chi_to_add += D[0]*dYdx*dYdx/((double) n_overlapped);
                        }
                    }
                }
                
                Chi_to_add = Chi_to_add*dx[0];
                Chi_integrated_local += Chi_to_add;
                
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
            &Chi_integrated_local,
            &Chi_integrated_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the enstrophy integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Chi_integrated_global;
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double Chi_integrated_local = 0.0;
        double Chi_integrated_global = 0.0;
        
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
                 * Register the patch, mass fraction, pressure and temperature in the flow model
                 * and compute the corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to mass fraction, pressure and temperature data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                boost::shared_ptr<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getGlobalCellData("PRESSURE");
                
                boost::shared_ptr<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getGlobalCellData("TEMPERATURE");
                
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_1_pressure = num_ghosts_pressure[1];
                const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                const int num_ghosts_1_temperature = num_ghosts_temperature[1];
                const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                
                double Chi_to_add = 0.0;
                
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
                    
                    if (false) // (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*4)
                    {
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
                                
                                // Compute indices of current and neighboring cells of mass fraction.
                                const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_LLLL = (relative_idx_lo_0 + i - 4 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_LLL = (relative_idx_lo_0 + i - 3 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_LL = (relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_RR = (relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_RRR = (relative_idx_lo_0 + i + 3 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_RRRR = (relative_idx_lo_0 + i + 4 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_BBBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 4 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_TT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_TTTT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 4 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                // Compute linear indices of pressure and temperature.
                                const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                
                                const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                
                                const double dYdx = (-1.0/280.0*Y[0][idx_x_RRRR] + 4.0/105.0*Y[0][idx_x_RRR] -
                                                     1.0/5.0*Y[0][idx_x_RR] + 4.0/5.0*Y[0][idx_x_R] -
                                                     4.0/5.0*Y[0][idx_x_L] + 1.0/5.0*Y[0][idx_x_LL] -
                                                     4.0/105.0*Y[0][idx_x_LLL] + 1.0/280.0*Y[0][idx_x_LLLL])/dx[0];
                                
                                const double dYdy = (-1.0/280.0*Y[0][idx_y_TTTT] + 4.0/105.0*Y[0][idx_y_TTT] -
                                                     1.0/5.0*Y[0][idx_y_TT] + 4.0/5.0*Y[0][idx_y_T] -
                                                     4.0/5.0*Y[0][idx_y_B] + 1.0/5.0*Y[0][idx_y_BB] -
                                                     4.0/105.0*Y[0][idx_y_BBB] + 1.0/280.0*Y[0][idx_y_BBBB])/dx[1];
                                
                                std::vector<double> D;
                                D.resize(d_num_species);
                                
                                std::vector<double*> D_ptr;
                                D_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    D_ptr.push_back(&D[si]);
                                }
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                                }
                                
                                d_equation_of_mass_diffusivity_mixing_rules->
                                    getMassDiffusivities(
                                        D_ptr,
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                Chi_to_add += D[0]*(dYdx*dYdx + dYdy*dYdy)/((double) n_overlapped);
                            }
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*3)
                    {
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
                                
                                // Compute indices of current and neighboring cells of mass fraction.
                                const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_LLL = (relative_idx_lo_0 + i - 3 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_LL = (relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_RR = (relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_RRR = (relative_idx_lo_0 + i + 3 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_TT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                // Compute linear indices of pressure and temperature.
                                const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                
                                const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                
                                const double dYdx = (1.0/60.0*Y[0][idx_x_RRR] - 3.0/20.0*Y[0][idx_x_RR] +
                                                     3.0/4.0*Y[0][idx_x_R] - 3.0/4.0*Y[0][idx_x_L] +
                                                     3.0/20.0*Y[0][idx_x_LL] - 1.0/60.0*Y[0][idx_x_LLL])/dx[0];
                                
                                const double dYdy = (1.0/60.0*Y[0][idx_y_TTT] - 3.0/20.0*Y[0][idx_y_TT] +
                                                     3.0/4.0*Y[0][idx_y_T] - 3.0/4.0*Y[0][idx_y_B] +
                                                     3.0/20.0*Y[0][idx_y_BB] - 1.0/60.0*Y[0][idx_y_BBB])/dx[1];
                                
                                std::vector<double> D;
                                D.resize(d_num_species);
                                
                                std::vector<double*> D_ptr;
                                D_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    D_ptr.push_back(&D[si]);
                                }
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                                }
                                
                                d_equation_of_mass_diffusivity_mixing_rules->
                                    getMassDiffusivities(
                                        D_ptr,
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                Chi_to_add += D[0]*(dYdx*dYdx + dYdy*dYdy)/((double) n_overlapped);
                            }
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*2)
                    {
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
                                
                                // Compute indices of current and neighboring cells of mass fraction.
                                const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_LL = (relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_RR = (relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_TT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                // Compute linear indices of pressure and temperature.
                                const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                
                                const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                
                                const double dYdx = (-1.0/12.0*Y[0][idx_x_RR] + 2.0/3.0*Y[0][idx_x_R] -
                                                     2.0/3.0*Y[0][idx_x_L] + 1.0/12.0*Y[0][idx_x_LL])/dx[0];
                                
                                const double dYdy = (-1.0/12.0*Y[0][idx_y_TT] + 2.0/3.0*Y[0][idx_y_T] -
                                                     2.0/3.0*Y[0][idx_y_B] + 1.0/12.0*Y[0][idx_y_BB])/dx[1];
                                
                                std::vector<double> D;
                                D.resize(d_num_species);
                                
                                std::vector<double*> D_ptr;
                                D_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    D_ptr.push_back(&D[si]);
                                }
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                                }
                                
                                d_equation_of_mass_diffusivity_mixing_rules->
                                    getMassDiffusivities(
                                        D_ptr,
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                Chi_to_add += D[0]*(dYdx*dYdx + dYdy*dYdy)/((double) n_overlapped);
                            }
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim))
                    {
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
                                
                                // Compute indices of current and neighboring cells of mass fraction.
                                const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                // Compute linear indices of pressure and temperature.
                                const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                
                                const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                
                                const double dYdx = (0.5*Y[0][idx_x_R] - 0.5*Y[0][idx_x_L])/dx[0];
                                const double dYdy = (0.5*Y[0][idx_y_T] - 0.5*Y[0][idx_y_B])/dx[1];
                                
                                std::vector<double> D;
                                D.resize(d_num_species);
                                
                                std::vector<double*> D_ptr;
                                D_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    D_ptr.push_back(&D[si]);
                                }
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                                }
                                
                                d_equation_of_mass_diffusivity_mixing_rules->
                                    getMassDiffusivities(
                                        D_ptr,
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                Chi_to_add += D[0]*(dYdx*dYdx + dYdy*dYdy)/((double) n_overlapped);
                            }
                        }
                    }
                    else
                    {
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
                                
                                // Compute indices of current and neighboring cells of mass fraction.
                                const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                // Compute linear indices of pressure and temperature.
                                const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                                
                                const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                                
                                double dYdx, dYdy;
                                
                                if (i == -num_ghosts_0_mass_fraction)
                                {
                                    dYdx = (Y[0][idx_x_R] - Y[0][idx_mass_fraction])/dx[0];
                                }
                                else if (i == interior_dim_0 + num_ghosts_0_mass_fraction - 1)
                                {
                                    dYdx = (Y[0][idx_mass_fraction] - Y[0][idx_x_L])/dx[0];
                                }
                                else
                                {
                                    dYdx = (0.5*Y[0][idx_x_R] - 0.5*Y[0][idx_x_L])/dx[0];
                                }
                                
                                if (j == -num_ghosts_1_mass_fraction)
                                {
                                    dYdy = (Y[0][idx_y_T] - Y[0][idx_mass_fraction])/dx[1];
                                }
                                else if (j == interior_dim_1 + num_ghosts_1_mass_fraction - 1)
                                {
                                    dYdy = (Y[0][idx_mass_fraction] - Y[0][idx_y_B])/dx[1];
                                }
                                else
                                {
                                    dYdy = (0.5*Y[0][idx_y_T] - 0.5*Y[0][idx_y_B])/dx[1];
                                }
                                
                                std::vector<double> D;
                                D.resize(d_num_species);
                                
                                std::vector<double*> D_ptr;
                                D_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    D_ptr.push_back(&D[si]);
                                }
                                
                                std::vector<const double*> Y_ptr;
                                Y_ptr.reserve(d_num_species);
                                for (int si = 0; si < d_num_species; si++)
                                {
                                    Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                                }
                                
                                d_equation_of_mass_diffusivity_mixing_rules->
                                    getMassDiffusivities(
                                        D_ptr,
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                Chi_to_add += D[0]*(dYdx*dYdx + dYdy*dYdy)/((double) n_overlapped);
                            }
                        }
                    }
                }
                
                Chi_to_add = Chi_to_add*dx[0]*dx[1];
                Chi_integrated_local += Chi_to_add;
                
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
            &Chi_integrated_local,
            &Chi_integrated_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the enstrophy integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Chi_integrated_global;
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double Chi_integrated_local = 0.0;
        double Chi_integrated_global = 0.0;
        
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
                 * Register the patch, mass fraction, pressure and temperature in the flow model
                 * and compute the corresponding cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to mass fraction, pressure and temperature data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                boost::shared_ptr<pdat::CellData<double> > data_pressure =
                    d_flow_model_tmp->getGlobalCellData("PRESSURE");
                
                boost::shared_ptr<pdat::CellData<double> > data_temperature =
                    d_flow_model_tmp->getGlobalCellData("TEMPERATURE");
                
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                double* p = data_pressure->getPointer(0);
                double* T = data_temperature->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int num_ghosts_2_mass_fraction = num_ghosts_mass_fraction[2];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                const int ghostcell_dim_1_mass_fraction = ghostcell_dims_mass_fraction[1];
                
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_1_pressure = num_ghosts_pressure[1];
                const int num_ghosts_2_pressure = num_ghosts_pressure[2];
                const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                const int ghostcell_dim_1_pressure = ghostcell_dims_pressure[1];
                
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                const int num_ghosts_1_temperature = num_ghosts_temperature[1];
                const int num_ghosts_2_temperature = num_ghosts_temperature[2];
                const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                const int ghostcell_dim_1_temperature = ghostcell_dims_temperature[1];
                
                double Chi_to_add = 0.0;
                
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
                    
                    if (false) // (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*4)
                    {
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
                                    
                                    // Compute indices of current and neighboring cells.
                                    const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_LLLL = (relative_idx_lo_0 + i - 4 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_LLL = (relative_idx_lo_0 + i - 3 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_LL = (relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_RR = (relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_RRR = (relative_idx_lo_0 + i + 3 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_RRRR = (relative_idx_lo_0 + i + 4 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_BBBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 4 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_TT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_TTTT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 4 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_BBBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 4 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 3 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 2 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_F = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_FF = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 2 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 3 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_FFFF = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 4 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    // Compute linear indices of pressure and temperature.
                                    const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                        (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                            ghostcell_dim_1_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                        (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                            ghostcell_dim_1_temperature;
                                    
                                    const double dYdx = (-1.0/280.0*Y[0][idx_x_RRRR] + 4.0/105.0*Y[0][idx_x_RRR] -
                                                         1.0/5.0*Y[0][idx_x_RR] + 4.0/5.0*Y[0][idx_x_R] -
                                                         4.0/5.0*Y[0][idx_x_L] + 1.0/5.0*Y[0][idx_x_LL] -
                                                         4.0/105.0*Y[0][idx_x_LLL] + 1.0/280.0*Y[0][idx_x_LLLL])/dx[0];
                                    
                                    const double dYdy = (-1.0/280.0*Y[0][idx_y_TTTT] + 4.0/105.0*Y[0][idx_y_TTT] -
                                                         1.0/5.0*Y[0][idx_y_TT] + 4.0/5.0*Y[0][idx_y_T] -
                                                         4.0/5.0*Y[0][idx_y_B] + 1.0/5.0*Y[0][idx_y_BB] -
                                                         4.0/105.0*Y[0][idx_y_BBB] + 1.0/280.0*Y[0][idx_y_BBBB])/dx[1];
                                    
                                    const double dYdz = (-1.0/280.0*Y[0][idx_z_FFFF] + 4.0/105.0*Y[0][idx_z_FFF] -
                                                         1.0/5.0*Y[0][idx_z_FF] + 4.0/5.0*Y[0][idx_z_F] -
                                                         4.0/5.0*Y[0][idx_z_B] + 1.0/5.0*Y[0][idx_z_BB] -
                                                         4.0/105.0*Y[0][idx_z_BBB] + 1.0/280.0*Y[0][idx_z_BBBB])/dx[2];
                                    
                                    std::vector<double> D;
                                    D.resize(d_num_species);
                                    
                                    std::vector<double*> D_ptr;
                                    D_ptr.reserve(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        D_ptr.push_back(&D[si]);
                                    }
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.reserve(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                                    }
                                    
                                    d_equation_of_mass_diffusivity_mixing_rules->
                                        getMassDiffusivities(
                                            D_ptr,
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    Chi_to_add += D[0]*(dYdx*dYdx + dYdy*dYdy + dYdz*dYdz)/((double) n_overlapped);
                                }
                            }
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*3)
                    {
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
                                    
                                    // Compute indices of current and neighboring cells.
                                    const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_LLL = (relative_idx_lo_0 + i - 3 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_LL = (relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_RR = (relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_RRR = (relative_idx_lo_0 + i + 3 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_TT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 3 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 2 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_F = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_FF = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 2 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 3 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    // Compute linear indices of pressure and temperature.
                                    const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                        (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                            ghostcell_dim_1_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                        (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                            ghostcell_dim_1_temperature;
                                    
                                    const double dYdx = (1.0/60.0*Y[0][idx_x_RRR] - 3.0/20.0*Y[0][idx_x_RR] +
                                                         3.0/4.0*Y[0][idx_x_R] - 3.0/4.0*Y[0][idx_x_L] +
                                                         3.0/20.0*Y[0][idx_x_LL] - 1.0/60.0*Y[0][idx_x_LLL])/dx[0];
                                    
                                    const double dYdy = (1.0/60.0*Y[0][idx_y_TTT] - 3.0/20.0*Y[0][idx_y_TT] +
                                                         3.0/4.0*Y[0][idx_y_T] - 3.0/4.0*Y[0][idx_y_B] +
                                                         3.0/20.0*Y[0][idx_y_BB] - 1.0/60.0*Y[0][idx_y_BBB])/dx[1];
                                    
                                    const double dYdz = (1.0/60.0*Y[0][idx_z_FFF] - 3.0/20.0*Y[0][idx_z_FF] +
                                                         3.0/4.0*Y[0][idx_z_F] - 3.0/4.0*Y[0][idx_z_B] +
                                                         3.0/20.0*Y[0][idx_z_BB] - 1.0/60.0*Y[0][idx_z_BBB])/dx[2];
                                    
                                    std::vector<double> D;
                                    D.resize(d_num_species);
                                    
                                    std::vector<double*> D_ptr;
                                    D_ptr.reserve(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        D_ptr.push_back(&D[si]);
                                    }
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.reserve(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                                    }
                                    
                                    d_equation_of_mass_diffusivity_mixing_rules->
                                        getMassDiffusivities(
                                            D_ptr,
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    Chi_to_add += D[0]*(dYdx*dYdx + dYdy*dYdy + dYdz*dYdz)/((double) n_overlapped);
                                }
                            }
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*2)
                    {
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
                                    
                                    // Compute indices of current and neighboring cells.
                                    const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_LL = (relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_RR = (relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_TT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 2 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_F = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_FF = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 2 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    // Compute linear indices of pressure and temperature.
                                    const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                        (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                            ghostcell_dim_1_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                        (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                            ghostcell_dim_1_temperature;
                                    
                                    const double dYdx = (-1.0/12.0*Y[0][idx_x_RR] + 2.0/3.0*Y[0][idx_x_R] -
                                                         2.0/3.0*Y[0][idx_x_L] + 1.0/12.0*Y[0][idx_x_LL])/dx[0];
                                    
                                    const double dYdy = (-1.0/12.0*Y[0][idx_y_TT] + 2.0/3.0*Y[0][idx_y_T] -
                                                         2.0/3.0*Y[0][idx_y_B] + 1.0/12.0*Y[0][idx_y_BB])/dx[1];
                                    
                                    const double dYdz = (-1.0/12.0*Y[0][idx_z_FF] + 2.0/3.0*Y[0][idx_z_F] -
                                                         2.0/3.0*Y[0][idx_z_B] + 1.0/12.0*Y[0][idx_z_BB])/dx[2];
                                    
                                    std::vector<double> D;
                                    D.resize(d_num_species);
                                    
                                    std::vector<double*> D_ptr;
                                    D_ptr.reserve(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        D_ptr.push_back(&D[si]);
                                    }
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.reserve(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                                    }
                                    
                                    d_equation_of_mass_diffusivity_mixing_rules->
                                        getMassDiffusivities(
                                            D_ptr,
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    Chi_to_add += D[0]*(dYdx*dYdx + dYdy*dYdy + dYdz*dYdz)/((double) n_overlapped);
                                }
                            }
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim))
                    {
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
                                    
                                    // Compute indices of current and neighboring cells.
                                    const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_F = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    // Compute linear indices of pressure and temperature.
                                    const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                        (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                            ghostcell_dim_1_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                        (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                            ghostcell_dim_1_temperature;
                                    
                                    const double dYdx = (0.5*Y[0][idx_x_R] - 0.5*Y[0][idx_x_L])/dx[0];
                                    const double dYdy = (0.5*Y[0][idx_y_T] - 0.5*Y[0][idx_y_B])/dx[1];
                                    const double dYdz = (0.5*Y[0][idx_z_F] - 0.5*Y[0][idx_z_B])/dx[2];
                                    
                                    std::vector<double> D;
                                    D.resize(d_num_species);
                                    
                                    std::vector<double*> D_ptr;
                                    D_ptr.reserve(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        D_ptr.push_back(&D[si]);
                                    }
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.reserve(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                                    }
                                    
                                    d_equation_of_mass_diffusivity_mixing_rules->
                                        getMassDiffusivities(
                                            D_ptr,
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    Chi_to_add += D[0]*(dYdx*dYdx + dYdy*dYdy + dYdz*dYdz)/((double) n_overlapped);
                                }
                            }
                        }
                    }
                    else
                    {
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
                                    
                                    // Compute indices of current and neighboring cells.
                                    const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_F = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    // Compute linear indices of pressure and temperature.
                                    const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                        (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                            ghostcell_dim_1_pressure;
                                    
                                    const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                        (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                            ghostcell_dim_1_temperature;
                                    
                                    double dYdx, dYdy, dYdz;
                                    
                                    if (i == -num_ghosts_0_mass_fraction)
                                    {
                                        dYdx = (Y[0][idx_x_R] - Y[0][idx_mass_fraction])/dx[0];
                                    }
                                    else if (i == interior_dim_0 + num_ghosts_0_mass_fraction - 1)
                                    {
                                        dYdx = (Y[0][idx_mass_fraction] - Y[0][idx_x_L])/dx[0];
                                    }
                                    else
                                    {
                                        dYdx = (0.5*Y[0][idx_x_R] - 0.5*Y[0][idx_x_L])/dx[0];
                                    }
                                    
                                    if (j == -num_ghosts_1_mass_fraction)
                                    {
                                        dYdy = (Y[0][idx_y_T] - Y[0][idx_mass_fraction])/dx[1];
                                    }
                                    else if (j == interior_dim_1 + num_ghosts_1_mass_fraction - 1)
                                    {
                                        dYdy = (Y[0][idx_mass_fraction] - Y[0][idx_y_B])/dx[1];
                                    }
                                    else
                                    {
                                        dYdy = (0.5*Y[0][idx_y_T] - 0.5*Y[0][idx_y_B])/dx[1];
                                    }
                                    
                                    if (k == -num_ghosts_2_mass_fraction)
                                    {
                                        dYdz = (Y[0][idx_z_F] - Y[0][idx_mass_fraction])/dx[2];
                                    }
                                    else if (k == interior_dim_2 + num_ghosts_2_mass_fraction - 1)
                                    {
                                        dYdz = (Y[0][idx_mass_fraction] - Y[0][idx_z_B])/dx[2];
                                    }
                                    else
                                    {
                                        dYdz = (0.5*Y[0][idx_z_F] - 0.5*Y[0][idx_z_B])/dx[2];
                                    }
                                    
                                    std::vector<double> D;
                                    D.resize(d_num_species);
                                    
                                    std::vector<double*> D_ptr;
                                    D_ptr.reserve(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        D_ptr.push_back(&D[si]);
                                    }
                                    
                                    std::vector<const double*> Y_ptr;
                                    Y_ptr.reserve(d_num_species);
                                    for (int si = 0; si < d_num_species; si++)
                                    {
                                        Y_ptr.push_back(&Y[si][idx_mass_fraction]);
                                    }
                                    
                                    d_equation_of_mass_diffusivity_mixing_rules->
                                        getMassDiffusivities(
                                            D_ptr,
                                            &p[idx_pressure],
                                            &T[idx_temperature],
                                            Y_ptr);
                                    
                                    Chi_to_add += D[0]*(dYdx*dYdx + dYdy*dYdy + dYdz*dYdz)/((double) n_overlapped);
                                }
                            }
                        }
                    }
                }
                
                Chi_to_add = Chi_to_add*dx[0]*dx[1]*dx[2];
                Chi_integrated_local += Chi_to_add;
                
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
            &Chi_integrated_local,
            &Chi_integrated_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the enstrophy integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Chi_integrated_global;
        }
    }
}


/*
 * Output numerical interface thickness to a file.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::outputNumericalInterfaceThickness(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
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
        f_out.open(stat_dump_filename.c_str(), std::ios::app);
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
    
    boost::shared_ptr<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
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
    
    if (d_dim == tbox::Dimension(1))
    {
        double grad_mag_max_local = 0.0;
        double grad_mag_max_global = 0.0;
        
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
                 * Register the patch, mass fraction in the flow model and compute the corresponding
                 * cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                
                for (hier::BoxContainer::BoxContainerConstIterator ib(patch_visible_boxes.begin());
                     ib != patch_visible_boxes.end();
                     ib++)
                {
                    const hier::Box& patch_visible_box = *ib;
                    
                    const hier::IntVector interior_dims = patch_visible_box.numberCells();
                    
                    const int interior_dim_0 = interior_dims[0];
                    
                    const hier::Index& index_lo = patch_visible_box.lower();
                    const hier::Index relative_index_lo = index_lo - patch_index_lo;
                    
                    const int relative_idx_lo_0 = relative_index_lo[0];
                    
                    if (false) // (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*4)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells of mass fraction.
                            const int idx_x_LLLL = relative_idx_lo_0 + i - 4 + num_ghosts_0_mass_fraction;
                            const int idx_x_LLL = relative_idx_lo_0 + i - 3 + num_ghosts_0_mass_fraction;
                            const int idx_x_LL = relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction;
                            const int idx_x_L = relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_R = relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_RR = relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction;
                            const int idx_x_RRR = relative_idx_lo_0 + i + 3 + num_ghosts_0_mass_fraction;
                            const int idx_x_RRRR = relative_idx_lo_0 + i + 4 + num_ghosts_0_mass_fraction;
                            
                            const double dYdx = (-1.0/280.0*Y[0][idx_x_RRRR] + 4.0/105.0*Y[0][idx_x_RRR] -
                                                 1.0/5.0*Y[0][idx_x_RR] + 4.0/5.0*Y[0][idx_x_R] -
                                                 4.0/5.0*Y[0][idx_x_L] + 1.0/5.0*Y[0][idx_x_LL] -
                                                 4.0/105.0*Y[0][idx_x_LLL] + 1.0/280.0*Y[0][idx_x_LLLL])/dx[0];
                            
                            const double grad_mag = fabs(dYdx);
                            grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*3)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells of mass fraction.
                            const int idx_x_LLL = relative_idx_lo_0 + i - 3 + num_ghosts_0_mass_fraction;
                            const int idx_x_LL = relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction;
                            const int idx_x_L = relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_R = relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_RR = relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction;
                            const int idx_x_RRR = relative_idx_lo_0 + i + 3 + num_ghosts_0_mass_fraction;
                            
                            const double dYdx = (1.0/60.0*Y[0][idx_x_RRR] - 3.0/20.0*Y[0][idx_x_RR] +
                                                 3.0/4.0*Y[0][idx_x_R] - 3.0/4.0*Y[0][idx_x_L] +
                                                 3.0/20.0*Y[0][idx_x_LL] - 1.0/60.0*Y[0][idx_x_LLL])/dx[0];
                            
                            const double grad_mag = fabs(dYdx);
                            grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*2)
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells of mass fraction.
                            const int idx_x_LL = relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction;
                            const int idx_x_L = relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_R = relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_RR = relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction;
                            
                            const double dYdx = (-1.0/12.0*Y[0][idx_x_RR] + 2.0/3.0*Y[0][idx_x_R] -
                                                 2.0/3.0*Y[0][idx_x_L] + 1.0/12.0*Y[0][idx_x_LL])/dx[0];
                            
                            const double grad_mag = fabs(dYdx);
                            grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim))
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells of mass fraction.
                            const int idx_x_L = relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_R = relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction;
                            
                            const double dYdx = (0.5*Y[0][idx_x_R] - 0.5*Y[0][idx_x_L])/dx[0];
                            
                            const double grad_mag = fabs(dYdx);
                            grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute indices of current and neighboring cells of mass fraction.
                            const int idx_mass_fraction = relative_idx_lo_0 + i + num_ghosts_0_mass_fraction;
                            const int idx_x_L = relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction;
                            const int idx_x_R = relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction;
                            
                            double dYdx;
                            
                            if (i == -num_ghosts_0_mass_fraction)
                            {
                                dYdx = (Y[0][idx_x_R] - Y[0][idx_mass_fraction])/dx[0];
                            }
                            else if (i == interior_dim_0 + num_ghosts_0_mass_fraction - 1)
                            {
                                dYdx = (Y[0][idx_mass_fraction] - Y[0][idx_x_L])/dx[0];
                            }
                            else
                            {
                                dYdx = (0.5*Y[0][idx_x_R] - 0.5*Y[0][idx_x_L])/dx[0];
                            }
                            
                            const double grad_mag = fabs(dYdx);
                            grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
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
         * Reduction to get the global integral.
         */
        
        mpi.Reduce(
            &grad_mag_max_local,
            &grad_mag_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0);
        
        /*
         * Output the numerical interface thickness (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            
            const double L_x = x_hi[0] - x_lo[0];
            
            const double dx_finest = L_x/finest_level_dims[0];
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << 1.0/(dx_finest*grad_mag_max_global);
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double grad_mag_max_local = 0.0;
        double grad_mag_max_global = 0.0;
        
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
                 * Register the patch, mass fraction in the flow model and compute the corresponding
                 * cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                
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
                    
                    if (false) // (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*4)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute indices of current and neighboring cells of mass fraction.
                                const int idx_x_LLLL = (relative_idx_lo_0 + i - 4 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_LLL = (relative_idx_lo_0 + i - 3 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_LL = (relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_RR = (relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_RRR = (relative_idx_lo_0 + i + 3 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_RRRR = (relative_idx_lo_0 + i + 4 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_BBBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 4 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_TT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_TTTT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 4 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const double dYdx = (-1.0/280.0*Y[0][idx_x_RRRR] + 4.0/105.0*Y[0][idx_x_RRR] -
                                                     1.0/5.0*Y[0][idx_x_RR] + 4.0/5.0*Y[0][idx_x_R] -
                                                     4.0/5.0*Y[0][idx_x_L] + 1.0/5.0*Y[0][idx_x_LL] -
                                                     4.0/105.0*Y[0][idx_x_LLL] + 1.0/280.0*Y[0][idx_x_LLLL])/dx[0];
                                
                                const double dYdy = (-1.0/280.0*Y[0][idx_y_TTTT] + 4.0/105.0*Y[0][idx_y_TTT] -
                                                     1.0/5.0*Y[0][idx_y_TT] + 4.0/5.0*Y[0][idx_y_T] -
                                                     4.0/5.0*Y[0][idx_y_B] + 1.0/5.0*Y[0][idx_y_BB] -
                                                     4.0/105.0*Y[0][idx_y_BBB] + 1.0/280.0*Y[0][idx_y_BBBB])/dx[1];
                                
                                const double grad_mag = sqrt(dYdx*dYdx + dYdy*dYdy);
                                grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
                            }
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*3)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute indices of current and neighboring cells of mass fraction.
                                const int idx_x_LLL = (relative_idx_lo_0 + i - 3 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_LL = (relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_RR = (relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_RRR = (relative_idx_lo_0 + i + 3 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_TT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const double dYdx = (1.0/60.0*Y[0][idx_x_RRR] - 3.0/20.0*Y[0][idx_x_RR] +
                                                     3.0/4.0*Y[0][idx_x_R] - 3.0/4.0*Y[0][idx_x_L] +
                                                     3.0/20.0*Y[0][idx_x_LL] - 1.0/60.0*Y[0][idx_x_LLL])/dx[0];
                                
                                const double dYdy = (1.0/60.0*Y[0][idx_y_TTT] - 3.0/20.0*Y[0][idx_y_TT] +
                                                     3.0/4.0*Y[0][idx_y_T] - 3.0/4.0*Y[0][idx_y_B] +
                                                     3.0/20.0*Y[0][idx_y_BB] - 1.0/60.0*Y[0][idx_y_BBB])/dx[1];
                                
                                const double grad_mag = sqrt(dYdx*dYdx + dYdy*dYdy);
                                grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
                            }
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*2)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute indices of current and neighboring cells of mass fraction.
                                const int idx_x_LL = (relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_RR = (relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_TT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const double dYdx = (-1.0/12.0*Y[0][idx_x_RR] + 2.0/3.0*Y[0][idx_x_R] -
                                                     2.0/3.0*Y[0][idx_x_L] + 1.0/12.0*Y[0][idx_x_LL])/dx[0];
                                
                                const double dYdy = (-1.0/12.0*Y[0][idx_y_TT] + 2.0/3.0*Y[0][idx_y_T] -
                                                     2.0/3.0*Y[0][idx_y_B] + 1.0/12.0*Y[0][idx_y_BB])/dx[1];
                                
                                const double grad_mag = sqrt(dYdx*dYdx + dYdy*dYdy);
                                grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
                            }
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim))
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute indices of current and neighboring cells of mass fraction.
                                const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const double dYdx = (0.5*Y[0][idx_x_R] - 0.5*Y[0][idx_x_L])/dx[0];
                                const double dYdy = (0.5*Y[0][idx_y_T] - 0.5*Y[0][idx_y_B])/dx[1];
                                
                                const double grad_mag = sqrt(dYdx*dYdx + dYdy*dYdy);
                                grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
                            }
                        }
                    }
                    else
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute indices of current and neighboring cells of mass fraction.
                                const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                    (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction;
                                
                                double dYdx, dYdy;
                                
                                if (i == -num_ghosts_0_mass_fraction)
                                {
                                    dYdx = (Y[0][idx_x_R] - Y[0][idx_mass_fraction])/dx[0];
                                }
                                else if (i == interior_dim_0 + num_ghosts_0_mass_fraction - 1)
                                {
                                    dYdx = (Y[0][idx_mass_fraction] - Y[0][idx_x_L])/dx[0];
                                }
                                else
                                {
                                    dYdx = (0.5*Y[0][idx_x_R] - 0.5*Y[0][idx_x_L])/dx[0];
                                }
                                
                                if (j == -num_ghosts_1_mass_fraction)
                                {
                                    dYdy = (Y[0][idx_y_T] - Y[0][idx_mass_fraction])/dx[1];
                                }
                                else if (j == interior_dim_1 + num_ghosts_1_mass_fraction - 1)
                                {
                                    dYdy = (Y[0][idx_mass_fraction] - Y[0][idx_y_B])/dx[1];
                                }
                                else
                                {
                                    dYdy = (0.5*Y[0][idx_y_T] - 0.5*Y[0][idx_y_B])/dx[1];
                                }
                                
                                const double grad_mag = sqrt(dYdx*dYdx + dYdy*dYdy);
                                grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
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
         * Reduction to get the global integral.
         */
        
        mpi.Reduce(
            &grad_mag_max_local,
            &grad_mag_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0);
        
        /*
         * Output the numerical interface thickness (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            
            const double L_x = x_hi[0] - x_lo[0];
            const double L_y = x_hi[1] - x_lo[1];
            
            double dx_finest = L_x/finest_level_dims[0];
            dx_finest = fmax(dx_finest, L_y/finest_level_dims[1]);
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << 1.0/(dx_finest*grad_mag_max_global);
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double grad_mag_max_local = 0.0;
        double grad_mag_max_global = 0.0;
        
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
                 * Register the patch, mass fraction in the flow model and compute the corresponding
                 * cell data.
                 */
                
                d_flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = d_flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTION", num_ghosts));
                
                d_flow_model_tmp->registerDerivedCellVariable(num_subghosts_of_data);
                
                d_flow_model_tmp->computeGlobalDerivedCellData();
                
                /*
                 * Get the pointers to mass fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mass_fraction =
                    d_flow_model_tmp->getGlobalCellData("MASS_FRACTION");
                
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fraction->getPointer(si));
                }
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fraction = data_mass_fraction->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fraction = data_mass_fraction->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fraction = num_ghosts_mass_fraction[0];
                const int num_ghosts_1_mass_fraction = num_ghosts_mass_fraction[1];
                const int num_ghosts_2_mass_fraction = num_ghosts_mass_fraction[2];
                const int ghostcell_dim_0_mass_fraction = ghostcell_dims_mass_fraction[0];
                const int ghostcell_dim_1_mass_fraction = ghostcell_dims_mass_fraction[1];
                
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
                    
                    if (false) // (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*4)
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute indices of current and neighboring cells.
                                    const int idx_x_LLLL = (relative_idx_lo_0 + i - 4 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_LLL = (relative_idx_lo_0 + i - 3 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_LL = (relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_RR = (relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_RRR = (relative_idx_lo_0 + i + 3 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_RRRR = (relative_idx_lo_0 + i + 4 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_BBBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 4 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_TT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_TTTT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 4 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_BBBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 4 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 3 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 2 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_F = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_FF = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 2 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 3 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_FFFF = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 4 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                        
                                    const double dYdx = (-1.0/280.0*Y[0][idx_x_RRRR] + 4.0/105.0*Y[0][idx_x_RRR] -
                                                         1.0/5.0*Y[0][idx_x_RR] + 4.0/5.0*Y[0][idx_x_R] -
                                                         4.0/5.0*Y[0][idx_x_L] + 1.0/5.0*Y[0][idx_x_LL] -
                                                         4.0/105.0*Y[0][idx_x_LLL] + 1.0/280.0*Y[0][idx_x_LLLL])/dx[0];
                                    
                                    const double dYdy = (-1.0/280.0*Y[0][idx_y_TTTT] + 4.0/105.0*Y[0][idx_y_TTT] -
                                                         1.0/5.0*Y[0][idx_y_TT] + 4.0/5.0*Y[0][idx_y_T] -
                                                         4.0/5.0*Y[0][idx_y_B] + 1.0/5.0*Y[0][idx_y_BB] -
                                                         4.0/105.0*Y[0][idx_y_BBB] + 1.0/280.0*Y[0][idx_y_BBBB])/dx[1];
                                    
                                    const double dYdz = (-1.0/280.0*Y[0][idx_z_FFFF] + 4.0/105.0*Y[0][idx_z_FFF] -
                                                         1.0/5.0*Y[0][idx_z_FF] + 4.0/5.0*Y[0][idx_z_F] -
                                                         4.0/5.0*Y[0][idx_z_B] + 1.0/5.0*Y[0][idx_z_BB] -
                                                         4.0/105.0*Y[0][idx_z_BBB] + 1.0/280.0*Y[0][idx_z_BBBB])/dx[2];
                                    
                                    const double grad_mag = sqrt(dYdx*dYdx + dYdy*dYdy + dYdz*dYdz);
                                    grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
                                }
                            }
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*3)
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute indices of current and neighboring cells.
                                    const int idx_x_LLL = (relative_idx_lo_0 + i - 3 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_LL = (relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_RR = (relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_RRR = (relative_idx_lo_0 + i + 3 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_BBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_TT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_TTT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 3 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_BBB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 3 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 2 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_F = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_FF = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 2 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_FFF = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 3 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const double dYdx = (1.0/60.0*Y[0][idx_x_RRR] - 3.0/20.0*Y[0][idx_x_RR] +
                                                         3.0/4.0*Y[0][idx_x_R] - 3.0/4.0*Y[0][idx_x_L] +
                                                         3.0/20.0*Y[0][idx_x_LL] - 1.0/60.0*Y[0][idx_x_LLL])/dx[0];
                                    
                                    const double dYdy = (1.0/60.0*Y[0][idx_y_TTT] - 3.0/20.0*Y[0][idx_y_TT] +
                                                         3.0/4.0*Y[0][idx_y_T] - 3.0/4.0*Y[0][idx_y_B] +
                                                         3.0/20.0*Y[0][idx_y_BB] - 1.0/60.0*Y[0][idx_y_BBB])/dx[1];
                                    
                                    const double dYdz = (1.0/60.0*Y[0][idx_z_FFF] - 3.0/20.0*Y[0][idx_z_FF] +
                                                         3.0/4.0*Y[0][idx_z_F] - 3.0/4.0*Y[0][idx_z_B] +
                                                         3.0/20.0*Y[0][idx_z_BB] - 1.0/60.0*Y[0][idx_z_BBB])/dx[2];
                                    
                                    const double grad_mag = sqrt(dYdx*dYdx + dYdy*dYdy + dYdz*dYdz);
                                    grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
                                }
                            }
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim)*2)
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute indices of current and neighboring cells.
                                    const int idx_x_LL = (relative_idx_lo_0 + i - 2 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_RR = (relative_idx_lo_0 + i + 2 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_TT = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 2 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_BB = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 2 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_F = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_FF = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 2 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const double dYdx = (-1.0/12.0*Y[0][idx_x_RR] + 2.0/3.0*Y[0][idx_x_R] -
                                                         2.0/3.0*Y[0][idx_x_L] + 1.0/12.0*Y[0][idx_x_LL])/dx[0];
                                    
                                    const double dYdy = (-1.0/12.0*Y[0][idx_y_TT] + 2.0/3.0*Y[0][idx_y_T] -
                                                         2.0/3.0*Y[0][idx_y_B] + 1.0/12.0*Y[0][idx_y_BB])/dx[1];
                                    
                                    const double dYdz = (-1.0/12.0*Y[0][idx_z_FF] + 2.0/3.0*Y[0][idx_z_F] -
                                                         2.0/3.0*Y[0][idx_z_B] + 1.0/12.0*Y[0][idx_z_BB])/dx[1];
                                    
                                    const double grad_mag = sqrt(dYdx*dYdx + dYdy*dYdy + dYdz*dYdz);
                                    grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
                                }
                            }
                        }
                    }
                    else if (num_ghosts_mass_fraction >= hier::IntVector::getOne(d_dim))
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute indices of current and neighboring cells.
                                    const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_F = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const double dYdx = (0.5*Y[0][idx_x_R] - 0.5*Y[0][idx_x_L])/dx[0];
                                    const double dYdy = (0.5*Y[0][idx_y_T] - 0.5*Y[0][idx_y_B])/dx[1];
                                    const double dYdz = (0.5*Y[0][idx_z_F] - 0.5*Y[0][idx_z_B])/dx[2];
                                    
                                    const double grad_mag = sqrt(dYdx*dYdx + dYdy*dYdy + dYdz*dYdz);
                                    grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int k = 0; k < interior_dim_2; k++)
                        {
                            for (int j = 0; j < interior_dim_1; j++)
                            {
                                for (int i = 0; i < interior_dim_0; i++)
                                {
                                    // Compute indices of current and neighboring cells.
                                    const int idx_mass_fraction = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_L = (relative_idx_lo_0 + i - 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_x_R = (relative_idx_lo_0 + i + 1 + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j - 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_y_T = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + 1 + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_B = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k - 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    const int idx_z_F = (relative_idx_lo_0 + i + num_ghosts_0_mass_fraction) +
                                        (relative_idx_lo_1 + j + num_ghosts_1_mass_fraction)*ghostcell_dim_0_mass_fraction +
                                        (relative_idx_lo_2 + k + 1 + num_ghosts_2_mass_fraction)*ghostcell_dim_0_mass_fraction*
                                            ghostcell_dim_1_mass_fraction;
                                    
                                    double dYdx, dYdy, dYdz;
                                    
                                    if (i == -num_ghosts_0_mass_fraction)
                                    {
                                        dYdx = (Y[0][idx_x_R] - Y[0][idx_mass_fraction])/dx[0];
                                    }
                                    else if (i == interior_dim_0 + num_ghosts_0_mass_fraction - 1)
                                    {
                                        dYdx = (Y[0][idx_mass_fraction] - Y[0][idx_x_L])/dx[0];
                                    }
                                    else
                                    {
                                        dYdx = (0.5*Y[0][idx_x_R] - 0.5*Y[0][idx_x_L])/dx[0];
                                    }
                                    
                                    if (j == -num_ghosts_1_mass_fraction)
                                    {
                                        dYdy = (Y[0][idx_y_T] - Y[0][idx_mass_fraction])/dx[1];
                                    }
                                    else if (j == interior_dim_1 + num_ghosts_1_mass_fraction - 1)
                                    {
                                        dYdy = (Y[0][idx_mass_fraction] - Y[0][idx_y_B])/dx[1];
                                    }
                                    else
                                    {
                                        dYdy = (0.5*Y[0][idx_y_T] - 0.5*Y[0][idx_y_B])/dx[1];
                                    }
                                    
                                    if (k == -num_ghosts_2_mass_fraction)
                                    {
                                        dYdz = (Y[0][idx_z_F] - Y[0][idx_mass_fraction])/dx[2];
                                    }
                                    else if (k == interior_dim_2 + num_ghosts_2_mass_fraction - 1)
                                    {
                                        dYdz = (Y[0][idx_mass_fraction] - Y[0][idx_z_B])/dx[2];
                                    }
                                    else
                                    {
                                        dYdz = (0.5*Y[0][idx_z_F] - 0.5*Y[0][idx_z_B])/dx[2];
                                    }
                                    
                                    const double grad_mag = sqrt(dYdx*dYdx + dYdy*dYdy + dYdz*dYdz);
                                    grad_mag_max_local = fmax(grad_mag_max_local, grad_mag);
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
         * Reduction to get the global integral.
         */
        
        mpi.Reduce(
            &grad_mag_max_local,
            &grad_mag_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0);
        
        /*
         * Output the numerical interface thickness (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            const double* x_lo = d_grid_geometry->getXLower();
            const double* x_hi = d_grid_geometry->getXUpper();
            
            const double L_x = x_hi[0] - x_lo[0];
            const double L_y = x_hi[1] - x_lo[1];
            const double L_z = x_hi[2] - x_lo[2];
            
            double dx_finest = L_x/finest_level_dims[0];
            dx_finest = fmax(dx_finest, L_y/finest_level_dims[1]);
            dx_finest = fmax(dx_finest, L_z/finest_level_dims[2]);
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << 1.0/(dx_finest*grad_mag_max_global);
        }
    }
}
