#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"

#include <fstream>


class SBIStatisticsUtilities
{
    public:
        SBIStatisticsUtilities(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_WEAK_PTR<FlowModel> flow_model,
            const HAMERS_SHARED_PTR<EquationOfStateMixingRules> equation_of_state_mixing_rules,
            const HAMERS_SHARED_PTR<EquationOfMassDiffusivityMixingRules> equation_of_mass_diffusivity_mixing_rules,
            const HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
            const HAMERS_SHARED_PTR<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
            const HAMERS_SHARED_PTR<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_species(num_species),
                d_flow_model(flow_model),
                d_equation_of_state_mixing_rules(equation_of_state_mixing_rules),
                d_equation_of_mass_diffusivity_mixing_rules(equation_of_mass_diffusivity_mixing_rules),
                d_equation_of_shear_viscosity_mixing_rules(equation_of_shear_viscosity_mixing_rules),
                d_equation_of_bulk_viscosity_mixing_rules(equation_of_bulk_viscosity_mixing_rules),
                d_equation_of_thermal_conductivity_mixing_rules(equation_of_thermal_conductivity_mixing_rules),
                d_num_ghosts_derivative(3)
        {}
        
        /*
         * Output mixing width in x-direction to a file.
         */
        void
        outputMixingWidthInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output centroid in x-directioin to a file.
         */
        void
        outputCentroidInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output minimum interface location in x-direction to a file.
         */
        void
        outputInterfaceMinInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output maximum interface location in x-direction to a file.
         */
        void
        outputInterfaceMaxInXDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output minimum interface location in y-direction to a file.
         */
        void
        outputInterfaceMinInYDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output maximum interface location in y-direction to a file.
         */
        void
        outputInterfaceMaxInYDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output minimum interface location in z-direction to a file.
         */
        void
        outputInterfaceMinInZDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output maximum interface location in z-direction to a file.
         */
        void
        outputInterfaceMaxInZDirection(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output circulation to a file.
         */
        void
        outputCirculation(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output scalar dissipation rate of first species integrated to a file.
         */
        void
        outputScalarDissipationRateIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output enstrophy integrated to a file.
         */
        void
        outputEnstrophyIntegrated(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output number of cells to a file.
         */
        void
        outputNumberOfCells(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output weighted number of cells to a file.
         */
        void
        outputWeightedNumberOfCells(
            const std::string& stat_dump_filename,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * HAMERS_WEAK_PTR to FlowModel.
         */
        const HAMERS_WEAK_PTR<FlowModel> d_flow_model;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfStateMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfStateMixingRules>
            d_equation_of_state_mixing_rules;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfMassDiffusivityMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfMassDiffusivityMixingRules>
            d_equation_of_mass_diffusivity_mixing_rules;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfShearViscosityMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfShearViscosityMixingRules>
            d_equation_of_shear_viscosity_mixing_rules;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfBulkViscosityMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfBulkViscosityMixingRules>
            d_equation_of_bulk_viscosity_mixing_rules;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfThermalConductivityMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfThermalConductivityMixingRules>
            d_equation_of_thermal_conductivity_mixing_rules;
        
        /*
         * Number of ghost cells to use in taking derivatives.
         */
        const int d_num_ghosts_derivative;
        
};


/*
 * Output mixing width in x-direction to a file.
 */
void
SBIStatisticsUtilities::outputMixingWidthInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
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
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                
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
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_mass_fractions;
                        
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
                
                flow_model_tmp->unregisterPatch();
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
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
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
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                
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
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                            
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
                
                flow_model_tmp->unregisterPatch();
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
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
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
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index and grid spacing.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
                
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
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                        ghostcell_dim_1_mass_fractions;
                                
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
                
                flow_model_tmp->unregisterPatch();
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
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output centroid in x-direction to a file.
 */
void
SBIStatisticsUtilities::outputCentroidInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'CENTROID_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        double num_local = 0.0;
        double num_global = 0.0;
        
        double den_local = 0.0;
        double den_global = 0.0;
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_0 = dx[0];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_0 = x_lo_patch[0];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                
                double num_to_add = 0.0;
                double den_to_add = 0.0;
                
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
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_mass_fractions;
                        
                        const double x = (relative_idx_lo_0 + i + 0.5)*dx_0 + x_lo_patch_0;
                        
                        const double Y_corrected = fmax(fmin(Y[idx], 1.0), 0.0);
                        num_to_add += x*Y_corrected/((double) n_overlapped);
                        den_to_add += Y_corrected/((double) n_overlapped);
                    }
                }
                
                num_to_add = num_to_add*dx[0];
                den_to_add = den_to_add*dx[0];
                
                num_local += num_to_add;
                den_local += den_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global integrals.
         */
        
        mpi.Reduce(
            &num_local,
            &num_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            &den_local,
            &den_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the centroid (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            const double x_c = num_global/den_global;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << x_c;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double num_local = 0.0;
        double num_global = 0.0;
        
        double den_local = 0.0;
        double den_global = 0.0;
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_0 = dx[0];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_0 = x_lo_patch[0];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                
                double num_to_add = 0.0;
                double den_to_add = 0.0;
                
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
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                            
                            const double x = (relative_idx_lo_0 + i + 0.5)*dx_0 + x_lo_patch_0;
                            
                            const double Y_corrected = fmax(fmin(Y[idx], 1.0), 0.0);
                            num_to_add += x*Y_corrected/((double) n_overlapped);
                            den_to_add += Y_corrected/((double) n_overlapped);
                        }
                    }
                }
                
                num_to_add = num_to_add*dx[0]*dx[1];
                den_to_add = den_to_add*dx[0]*dx[1];
                
                num_local += num_to_add;
                den_local += den_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global integrals.
         */
        
        mpi.Reduce(
            &num_local,
            &num_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            &den_local,
            &den_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the centroid (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            const double x_c = num_global/den_global;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << x_c;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double num_local = 0.0;
        double num_global = 0.0;
        
        double den_local = 0.0;
        double den_global = 0.0;
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_0 = dx[0];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_0 = x_lo_patch[0];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
                
                double num_to_add = 0.0;
                double den_to_add = 0.0;
                
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
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                        ghostcell_dim_1_mass_fractions;
                                
                                const double x = (relative_idx_lo_0 + i + 0.5)*dx_0 + x_lo_patch_0;
                                
                                const double Y_corrected = fmax(fmin(Y[idx], 1.0), 0.0);
                                num_to_add += x*Y_corrected/((double) n_overlapped);
                                den_to_add += Y_corrected/((double) n_overlapped);
                            }
                        }
                    }
                }
                
                num_to_add = num_to_add*dx[0]*dx[1]*dx[2];
                den_to_add = den_to_add*dx[0]*dx[1]*dx[2];
                
                num_local += num_to_add;
                den_local += den_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global integrals.
         */
        
        mpi.Reduce(
            &num_local,
            &num_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            &den_local,
            &den_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the centroid (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            const double x_c = num_global/den_global;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << x_c;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
}


/*
 * Output minimum interface location in x-direction to a file.
 */
void
SBIStatisticsUtilities::outputInterfaceMinInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MIN_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the upper index of the physical domain.
     */
    
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        double interface_x_min_local = x_hi[0];
        double interface_x_min_global = x_hi[0];
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_0 = dx[0];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_0 = x_lo_patch[0];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                
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
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the linear index and the data to add.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_mass_fractions;
                        
                        const double x = (relative_idx_lo_0 + i + 0.5)*dx_0 + x_lo_patch_0;
                        
                        if (Y[idx] > 0.01 && Y[idx] < 0.99)
                        {
                            if (x < interface_x_min_local)
                            {
                                interface_x_min_local = x;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the minimum value.
         */
        
        mpi.Reduce(
            &interface_x_min_local,
            &interface_x_min_global,
            1,
            MPI_DOUBLE,
            MPI_MIN,
            0);
        
        /*
         * Compute and output the value (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << interface_x_min_global;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double interface_x_min_local = x_hi[0];
        double interface_x_min_global = x_hi[0];
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_0 = dx[0];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_0 = x_lo_patch[0];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                
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
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                            
                            const double x = (relative_idx_lo_0 + i + 0.5)*dx_0 + x_lo_patch_0;
                            
                            if (Y[idx] > 0.01 && Y[idx] < 0.99)
                            {
                                if (x < interface_x_min_local)
                                {
                                    interface_x_min_local = x;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the minimum value.
         */
        
        mpi.Reduce(
            &interface_x_min_local,
            &interface_x_min_global,
            1,
            MPI_DOUBLE,
            MPI_MIN,
            0);
        
        /*
         * Compute and output the value (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << interface_x_min_global;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double interface_x_min_local = x_hi[0];
        double interface_x_min_global = x_hi[0];
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_0 = dx[0];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_0 = x_lo_patch[0];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
                
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
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                        ghostcell_dim_1_mass_fractions;
                                
                                const double x = (relative_idx_lo_0 + i + 0.5)*dx_0 + x_lo_patch_0;
                                
                                if (Y[idx] > 0.01 && Y[idx] < 0.99)
                                {
                                    if (x < interface_x_min_local)
                                    {
                                        interface_x_min_local = x;
                                    }
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the minimum value.
         */
        
        mpi.Reduce(
            &interface_x_min_local,
            &interface_x_min_global,
            1,
            MPI_DOUBLE,
            MPI_MIN,
            0);
        
        /*
         * Compute and output the value (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << interface_x_min_global;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
}


/*
 * Output maximum interface location in x-direction to a file.
 */
void
SBIStatisticsUtilities::outputInterfaceMaxInXDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MAX_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the lower index of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    
    if (d_dim == tbox::Dimension(1))
    {
        double interface_x_max_local = x_lo[0];
        double interface_x_max_global = x_lo[0];
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_0 = dx[0];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_0 = x_lo_patch[0];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                
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
                    
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        /*
                         * Compute the linear index and the data to add.
                         */
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_mass_fractions;
                        
                        const double x = (relative_idx_lo_0 + i + 0.5)*dx_0 + x_lo_patch_0;
                        
                        if (Y[idx] > 0.01 && Y[idx] < 0.99)
                        {
                            if (x > interface_x_max_local)
                            {
                                interface_x_max_local = x;
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the maximum value.
         */
        
        mpi.Reduce(
            &interface_x_max_local,
            &interface_x_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0);
        
        /*
         * Compute and output the value (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << interface_x_max_global;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double interface_x_max_local = x_lo[0];
        double interface_x_max_global = x_lo[0];
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_0 = dx[0];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_0 = x_lo_patch[0];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                
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
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                            
                            const double x = (relative_idx_lo_0 + i + 0.5)*dx_0 + x_lo_patch_0;
                            
                            if (Y[idx] > 0.01 && Y[idx] < 0.99)
                            {
                                if (x > interface_x_max_local)
                                {
                                    interface_x_max_local = x;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the maximum value.
         */
        
        mpi.Reduce(
            &interface_x_max_local,
            &interface_x_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0);
        
        /*
         * Compute and output the value (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << interface_x_max_global;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double interface_x_max_local = x_lo[0];
        double interface_x_max_global = x_lo[0];
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_0 = dx[0];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_0 = x_lo_patch[0];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
                
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
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                        ghostcell_dim_1_mass_fractions;
                                
                                const double x = (relative_idx_lo_0 + i + 0.5)*dx_0 + x_lo_patch_0;
                                
                                if (Y[idx] > 0.01 && Y[idx] < 0.99)
                                {
                                    if (x > interface_x_max_local)
                                    {
                                        interface_x_max_local = x;
                                    }
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the maximum value.
         */
        
        mpi.Reduce(
            &interface_x_max_local,
            &interface_x_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0);
        
        /*
         * Compute and output the value (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << interface_x_max_global;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
}


/*
 * Output minimum interface location in y-direction to a file.
 */
void
SBIStatisticsUtilities::outputInterfaceMinInYDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MIN_Y' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the upper index of the physical domain.
     */
    
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'INTERFACE_MIN_Y' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double interface_y_min_local = x_hi[1];
        double interface_y_min_global = x_hi[1];
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_1 = dx[1];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_1 = x_lo_patch[1];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                
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
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                            
                            const double y = (relative_idx_lo_1 + j + 0.5)*dx_1 + x_lo_patch_1;
                            
                            if (Y[idx] > 0.01 && Y[idx] < 0.99)
                            {
                                if (y < interface_y_min_local)
                                {
                                    interface_y_min_local = y;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the minimum value.
         */
        
        mpi.Reduce(
            &interface_y_min_local,
            &interface_y_min_global,
            1,
            MPI_DOUBLE,
            MPI_MIN,
            0);
        
        /*
         * Compute and output the value (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << interface_y_min_global;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double interface_y_min_local = x_hi[1];
        double interface_y_min_global = x_hi[1];
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_1 = dx[1];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_1 = x_lo_patch[1];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
                
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
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                        ghostcell_dim_1_mass_fractions;
                                
                                const double y = (relative_idx_lo_1 + j + 0.5)*dx_1 + x_lo_patch_1;
                                
                                if (Y[idx] > 0.01 && Y[idx] < 0.99)
                                {
                                    if (y < interface_y_min_local)
                                    {
                                        interface_y_min_local = y;
                                    }
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the minimum value.
         */
        
        mpi.Reduce(
            &interface_y_min_local,
            &interface_y_min_global,
            1,
            MPI_DOUBLE,
            MPI_MIN,
            0);
        
        /*
         * Compute and output the value (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << interface_y_min_global;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
}


/*
 * Output maximum interface location in y-direction to a file.
 */
void
SBIStatisticsUtilities::outputInterfaceMaxInYDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MAX_Y' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the lower index of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'INTERFACE_MIN_Y' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double interface_y_max_local = x_lo[1];
        double interface_y_max_global = x_lo[1];
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_1 = dx[1];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_1 = x_lo_patch[1];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                
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
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            /*
                             * Compute the linear index and the data to add.
                             */
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                            
                            const double y = (relative_idx_lo_1 + j + 0.5)*dx_1 + x_lo_patch_1;
                            
                            if (Y[idx] > 0.01 && Y[idx] < 0.99)
                            {
                                if (y > interface_y_max_local)
                                {
                                    interface_y_max_local = y;
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the maximum value.
         */
        
        mpi.Reduce(
            &interface_y_max_local,
            &interface_y_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0);
        
        /*
         * Compute and output the value (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << interface_y_max_global;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double interface_y_max_local = x_lo[1];
        double interface_y_max_global = x_lo[1];
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_1 = dx[1];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_1 = x_lo_patch[1];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
                
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
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                        ghostcell_dim_1_mass_fractions;
                                
                                const double y = (relative_idx_lo_1 + j + 0.5)*dx_1 + x_lo_patch_1;
                                
                                if (Y[idx] > 0.01 && Y[idx] < 0.99)
                                {
                                    if (y > interface_y_max_local)
                                    {
                                        interface_y_max_local = y;
                                    }
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the maximum value.
         */
        
        mpi.Reduce(
            &interface_y_max_local,
            &interface_y_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0);
        
        /*
         * Compute and output the value (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << interface_y_max_global;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
}


/*
 * Output minimum interface location in z-direction to a file.
 */
void
SBIStatisticsUtilities::outputInterfaceMinInZDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MIN_Z' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the upper index of the physical domain.
     */
    
    const double* x_hi = d_grid_geometry->getXUpper();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'INTERFACE_MIN_Z' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'INTERFACE_MIN_Z' for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double interface_z_min_local = x_hi[2];
        double interface_z_min_global = x_hi[2];
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_2 = dx[2];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_2 = x_lo_patch[2];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
                
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
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                        ghostcell_dim_1_mass_fractions;
                                
                                const double z = (relative_idx_lo_2 + k + 0.5)*dx_2 + x_lo_patch_2;
                                
                                if (Y[idx] > 0.01 && Y[idx] < 0.99)
                                {
                                    if (z < interface_z_min_local)
                                    {
                                        interface_z_min_local = z;
                                    }
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the minimum value.
         */
        
        mpi.Reduce(
            &interface_z_min_local,
            &interface_z_min_global,
            1,
            MPI_DOUBLE,
            MPI_MIN,
            0);
        
        /*
         * Compute and output the value (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << interface_z_min_global;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
}


/*
 * Output maximum interface location in z-direction to a file.
 */
void
SBIStatisticsUtilities::outputInterfaceMaxInZDirection(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'INTERFACE_MAX_Z' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    /*
     * Get the lower index of the physical domain.
     */
    
    const double* x_lo = d_grid_geometry->getXLower();
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'INTERFACE_MIN_Z' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'INTERFACE_MIN_Z' for two-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double interface_z_max_local = x_lo[2];
        double interface_z_max_global = x_lo[2];
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch lower index, grid spacing and the lower spatial coordinate.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                const double dx_2 = dx[2];
                
                const double* const x_lo_patch = patch_geom->getXLower();
                
                const double x_lo_patch_2 = x_lo_patch[2];
                
                /*
                 * Register the patch and mass fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mass fraction data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                double* Y = data_mass_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
                
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
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                /*
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                        ghostcell_dim_1_mass_fractions;
                                
                                const double z = (relative_idx_lo_2 + k + 0.5)*dx_2 + x_lo_patch_2;
                                
                                if (Y[idx] > 0.01 && Y[idx] < 0.99)
                                {
                                    if (z > interface_z_max_local)
                                    {
                                        interface_z_max_local = z;
                                    }
                                }
                            }
                        }
                    }
                }
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the maximum value.
         */
        
        mpi.Reduce(
            &interface_z_max_local,
            &interface_z_max_global,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0);
        
        /*
         * Compute and output the value (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << interface_z_max_global;
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
}


/*
 * Output circulation to a file.
 */
void
SBIStatisticsUtilities::outputCirculation(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
        new ExtendedFlattenedHierarchy(
            *patch_hierarchy,
            0,
            num_levels - 1));
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'CIRCULATION' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double Circulation_local = 0.0;
        double Circulation_global = 0.0;
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch geometry.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector patch_dims = patch_box.numberCells();
                
                const int patch_dim_0 = patch_dims[0];
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and velocity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to velocity data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_velocity >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Circulation_to_add = 0.0;
                
                /*
                 * Initialize cell data for velocity derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue()*d_dim.getValue() - d_dim.getValue(),
                        hier::IntVector::getZero(d_dim)));
                
                double* dudy = data_velocity_derivatives->getPointer(0);
                double* dvdx = data_velocity_derivatives->getPointer(1);
                
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
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                        new DerivativeFirstOrder(
                            "first order derivative in x-direction",
                            d_dim,
                            DIRECTION::X_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                        new DerivativeFirstOrder(
                            "first order derivative in y-direction",
                            d_dim,
                            DIRECTION::Y_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    // Compute dudy.
                    derivative_first_order_y->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[1],
                        patch_visible_box,
                        0,
                        0);
                    
                    // Compute dvdx.
                    derivative_first_order_x->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[0],
                        patch_visible_box,
                        1,
                        1);
                    
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
                            const int idx = (relative_idx_lo_0 + i) +
                                (relative_idx_lo_1 + j)*patch_dim_0;
                            
                            const double omega = dvdx[idx] - dudy[idx];
                            Circulation_to_add += omega/((double) n_overlapped);
                        }
                    }
                }
                
                Circulation_to_add = Circulation_to_add*dx[0]*dx[1];
                Circulation_local += Circulation_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global integral.
         */
        
        mpi.Reduce(
            &Circulation_local,
            &Circulation_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the circulation (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << fabs(Circulation_global);
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double Circulation_y_local = 0.0;
        double Circulation_y_global = 0.0;
        
        double Circulation_z_local = 0.0;
        double Circulation_z_global = 0.0;
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch geometry.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector patch_dims = patch_box.numberCells();
                
                const int patch_dim_0 = patch_dims[0];
                const int patch_dim_1 = patch_dims[1];
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch and velocity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to velocity data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity = flow_model_tmp->getCellData("VELOCITY");
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_velocity >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Circulation_y_to_add = 0.0;
                double Circulation_z_to_add = 0.0;
                
                /*
                 * Initialize cell data for velocity derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue()*d_dim.getValue() - d_dim.getValue(),
                        hier::IntVector::getZero(d_dim)));
                
                double* dudy = data_velocity_derivatives->getPointer(0);
                double* dudz = data_velocity_derivatives->getPointer(1);
                double* dvdx = data_velocity_derivatives->getPointer(2);
                // double* dvdz = data_velocity_derivatives->getPointer(3);
                double* dwdx = data_velocity_derivatives->getPointer(4);
                // double* dwdy = data_velocity_derivatives->getPointer(5);
                
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
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                        new DerivativeFirstOrder(
                            "first order derivative in x-direction",
                            d_dim,
                            DIRECTION::X_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                        new DerivativeFirstOrder(
                            "first order derivative in y-direction",
                            d_dim,
                            DIRECTION::Y_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_z(
                        new DerivativeFirstOrder(
                            "first order derivative in z-direction",
                            d_dim,
                            DIRECTION::Z_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    // Compute dudy.
                    derivative_first_order_y->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[1],
                        patch_visible_box,
                        0,
                        0);
                    
                    // Compute dudz.
                    derivative_first_order_z->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[2],
                        patch_visible_box,
                        1,
                        0);
                    
                    // Compute dvdx.
                    derivative_first_order_x->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[0],
                        patch_visible_box,
                        2,
                        1);
                    
                    // Compute dvdz.
                    derivative_first_order_z->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[2],
                        patch_visible_box,
                        3,
                        1);
                    
                    // Compute dwdx.
                    derivative_first_order_x->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[0],
                        patch_visible_box,
                        4,
                        2);
                    
                    // Compute dwdy.
                    derivative_first_order_y->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[1],
                        patch_visible_box,
                        5,
                        2);
                    
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
                                
                                // Compute the linear indices.
                                const int idx = (relative_idx_lo_0 + i) +
                                    (relative_idx_lo_1 + j)*patch_dim_0 +
                                    (relative_idx_lo_2 + k)*patch_dim_0*
                                        patch_dim_1;
                                
                                const double omega_y = dudz[idx] - dwdx[idx];
                                const double omega_z = dvdx[idx] - dudy[idx];
                                Circulation_y_to_add += omega_y/((double) n_overlapped);
                                Circulation_z_to_add += omega_z/((double) n_overlapped);
                            }
                        }
                    }
                }
                
                Circulation_y_to_add = Circulation_y_to_add*dx[0]*dx[1]*dx[2];
                Circulation_z_to_add = Circulation_z_to_add*dx[0]*dx[1]*dx[2];
                
                Circulation_y_local += Circulation_y_to_add;
                Circulation_z_local += Circulation_z_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        /*
         * Reduction to get the global integral.
         */
        
        mpi.Reduce(
            &Circulation_y_local,
            &Circulation_y_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            &Circulation_z_local,
            &Circulation_z_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the circulation (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << 0.5*(fabs(Circulation_y_global) + fabs(Circulation_z_global));
        }
        
        if (mpi.getRank() == 0)
        {
            f_out.close();
        }
    }
}


/*
 * Output scalar dissipation rate of first species integrated to a file.
 */
void
SBIStatisticsUtilities::outputScalarDissipationRateIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch geometry.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch, mass fractions, pressure and temperature in the flow model
                 * and compute the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mass fraction, pressure and temperature data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    flow_model_tmp->getCellData("TEMPERATURE");
                
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fractions->getPointer(si));
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
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_mass_fractions >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Chi_to_add = 0.0;
                
                /*
                 * Initialize cell data for mass fraction derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
                
                double* dYdx = data_mass_fraction_derivatives->getPointer(0);
                
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
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                        new DerivativeFirstOrder(
                            "first order derivative in x-direction",
                            d_dim,
                            DIRECTION::X_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    // Compute dYdx.
                    derivative_first_order_x->computeDerivative(
                        data_mass_fraction_derivatives,
                        data_mass_fractions,
                        dx[0],
                        patch_visible_box,
                        0,
                        0);
                    
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
                        
                        // Compute linear indices of derivatives, mass fractions, pressure and temperature.
                        const int idx = relative_idx_lo_0 + i;
                        const int idx_mass_fractions = relative_idx_lo_0 + i + num_ghosts_0_mass_fractions;
                        const int idx_pressure = relative_idx_lo_0 + i + num_ghosts_0_pressure;
                        const int idx_temperature = relative_idx_lo_0 + i + num_ghosts_0_temperature;
                        
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
                            Y_ptr.push_back(&Y[si][idx_mass_fractions]);
                        }
                        
                        d_equation_of_mass_diffusivity_mixing_rules->
                            getMassDiffusivities(
                                D_ptr,
                                &p[idx_pressure],
                                &T[idx_temperature],
                                Y_ptr);
                        
                        Chi_to_add += D[0]*dYdx[idx]*dYdx[idx]/((double) n_overlapped);
                    }
                }
                
                Chi_to_add = Chi_to_add*dx[0];
                Chi_integrated_local += Chi_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
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
         * Output the scalar dissipation integral (only done by process 0).
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
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch geometry.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector patch_dims = patch_box.numberCells();
                
                const int patch_dim_0 = patch_dims[0];
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch, mass fractions, pressure and temperature in the flow model
                 * and compute the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mass fraction, pressure and temperature data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    flow_model_tmp->getCellData("TEMPERATURE");
                
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fractions->getPointer(si));
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
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                
                const int num_ghosts_0_pressure = num_ghosts_pressure[0];
                const int num_ghosts_1_pressure = num_ghosts_pressure[1];
                const int ghostcell_dim_0_pressure = ghostcell_dims_pressure[0];
                
                const int num_ghosts_0_temperature = num_ghosts_temperature[0];
                const int num_ghosts_1_temperature = num_ghosts_temperature[1];
                const int ghostcell_dim_0_temperature = ghostcell_dims_temperature[0];
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_mass_fractions >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Chi_to_add = 0.0;
                
                /*
                 * Initialize cell data for mass fraction derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
                
                double* dYdx = data_mass_fraction_derivatives->getPointer(0);
                double* dYdy = data_mass_fraction_derivatives->getPointer(1);
                
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
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                        new DerivativeFirstOrder(
                            "first order derivative in x-direction",
                            d_dim,
                            DIRECTION::X_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                        new DerivativeFirstOrder(
                            "first order derivative in y-direction",
                            d_dim,
                            DIRECTION::Y_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    // Compute dYdx.
                    derivative_first_order_x->computeDerivative(
                        data_mass_fraction_derivatives,
                        data_mass_fractions,
                        dx[0],
                        patch_visible_box,
                        0,
                        0);
                    
                    // Compute dYdy.
                    derivative_first_order_y->computeDerivative(
                        data_mass_fraction_derivatives,
                        data_mass_fractions,
                        dx[1],
                        patch_visible_box,
                        1,
                        0);
                    
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
                            
                            // Compute linear indices of derivatives, mass fractions, pressure and temperature.
                            const int idx = (relative_idx_lo_0 + i) +
                                (relative_idx_lo_1 + j)*patch_dim_0;
                            
                            const int idx_mass_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions;
                            
                            const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure;
                            
                            const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature;
                            
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
                                Y_ptr.push_back(&Y[si][idx_mass_fractions]);
                            }
                            
                            d_equation_of_mass_diffusivity_mixing_rules->
                                getMassDiffusivities(
                                    D_ptr,
                                    &p[idx_pressure],
                                    &T[idx_temperature],
                                    Y_ptr);
                            
                            Chi_to_add += D[0]*(dYdx[idx]*dYdx[idx] + dYdy[idx]*dYdy[idx])/((double) n_overlapped);
                        }
                    }
                }
                
                Chi_to_add = Chi_to_add*dx[0]*dx[1];
                Chi_integrated_local += Chi_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
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
         * Output the scalar dissipation integral (only done by process 0).
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
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch geometry.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector patch_dims = patch_box.numberCells();
                
                const int patch_dim_0 = patch_dims[0];
                const int patch_dim_1 = patch_dims[1];
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch, mass fractions, pressure and temperature in the flow model
                 * and compute the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MASS_FRACTIONS", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("PRESSURE", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("TEMPERATURE", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mass fraction, pressure and temperature data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fractions =
                    flow_model_tmp->getCellData("MASS_FRACTIONS");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_pressure =
                    flow_model_tmp->getCellData("PRESSURE");
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_temperature =
                    flow_model_tmp->getCellData("TEMPERATURE");
                
                std::vector<double*> Y;
                Y.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    Y.push_back(data_mass_fractions->getPointer(si));
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
                
                const hier::IntVector num_ghosts_mass_fractions = data_mass_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mass_fractions = data_mass_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_pressure = data_pressure->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_pressure = data_pressure->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_temperature = data_temperature->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_temperature = data_temperature->getGhostBox().numberCells();
                
                const int num_ghosts_0_mass_fractions = num_ghosts_mass_fractions[0];
                const int num_ghosts_1_mass_fractions = num_ghosts_mass_fractions[1];
                const int num_ghosts_2_mass_fractions = num_ghosts_mass_fractions[2];
                const int ghostcell_dim_0_mass_fractions = ghostcell_dims_mass_fractions[0];
                const int ghostcell_dim_1_mass_fractions = ghostcell_dims_mass_fractions[1];
                
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
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_mass_fractions >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Chi_to_add = 0.0;
                
                /*
                 * Initialize cell data for mass fraction derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_mass_fraction_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
                
                double* dYdx = data_mass_fraction_derivatives->getPointer(0);
                double* dYdy = data_mass_fraction_derivatives->getPointer(1);
                double* dYdz = data_mass_fraction_derivatives->getPointer(2);
                
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
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                        new DerivativeFirstOrder(
                            "first order derivative in x-direction",
                            d_dim,
                            DIRECTION::X_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                        new DerivativeFirstOrder(
                            "first order derivative in y-direction",
                            d_dim,
                            DIRECTION::Y_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_z(
                        new DerivativeFirstOrder(
                            "first order derivative in z-direction",
                            d_dim,
                            DIRECTION::Z_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    // Compute dYdx.
                    derivative_first_order_x->computeDerivative(
                        data_mass_fraction_derivatives,
                        data_mass_fractions,
                        dx[0],
                        patch_visible_box,
                        0,
                        0);
                    
                    // Compute dYdy.
                    derivative_first_order_y->computeDerivative(
                        data_mass_fraction_derivatives,
                        data_mass_fractions,
                        dx[1],
                        patch_visible_box,
                        1,
                        0);
                    
                    // Compute dYdz.
                    derivative_first_order_z->computeDerivative(
                        data_mass_fraction_derivatives,
                        data_mass_fractions,
                        dx[2],
                        patch_visible_box,
                        2,
                        0);
                    
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
                                
                                // Compute linear indices of derivatives, mass fractions, pressure and temperature.
                                const int idx = (relative_idx_lo_0 + i) +
                                    (relative_idx_lo_1 + j)*patch_dim_0 +
                                    (relative_idx_lo_2 + k)*patch_dim_0*
                                        patch_dim_1;
                                
                                const int idx_mass_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mass_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mass_fractions)*ghostcell_dim_0_mass_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mass_fractions)*ghostcell_dim_0_mass_fractions*
                                        ghostcell_dim_1_mass_fractions;
                                
                                const int idx_pressure = (relative_idx_lo_0 + i + num_ghosts_0_pressure) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_pressure)*ghostcell_dim_0_pressure +
                                    (relative_idx_lo_2 + k + num_ghosts_2_pressure)*ghostcell_dim_0_pressure*
                                        ghostcell_dim_1_pressure;
                                
                                const int idx_temperature = (relative_idx_lo_0 + i + num_ghosts_0_temperature) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_temperature)*ghostcell_dim_0_temperature +
                                    (relative_idx_lo_2 + k + num_ghosts_2_temperature)*ghostcell_dim_0_temperature*
                                        ghostcell_dim_1_temperature;
                                
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
                                    Y_ptr.push_back(&Y[si][idx_mass_fractions]);
                                }
                                
                                d_equation_of_mass_diffusivity_mixing_rules->
                                    getMassDiffusivities(
                                        D_ptr,
                                        &p[idx_pressure],
                                        &T[idx_temperature],
                                        Y_ptr);
                                
                                Chi_to_add += D[0]*(dYdx[idx]*dYdx[idx] + dYdy[idx]*dYdy[idx] + dYdz[idx]*dYdz[idx])/
                                    ((double) n_overlapped);
                            }
                        }
                    }
                }
                
                Chi_to_add = Chi_to_add*dx[0]*dx[1]*dx[2];
                Chi_integrated_local += Chi_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
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
         * Output the scalar dissipation integral (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Chi_integrated_global;
        }
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output enstrophy integrated to a file.
 */
void
SBIStatisticsUtilities::outputEnstrophyIntegrated(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
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
    
    HAMERS_SHARED_PTR<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
    
    HAMERS_SHARED_PTR<ExtendedFlattenedHierarchy> flattened_hierarchy(
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
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch geometry.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector patch_dims = patch_box.numberCells();
                
                const int patch_dim_0 = patch_dims[0];
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch, density and velocity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_density = flow_model_tmp->getCellData("DENSITY");
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity = flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                
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
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_velocity >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Omega_to_add = 0.0;
                
                /*
                 * Initialize cell data for velocity derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue()*d_dim.getValue() - d_dim.getValue(),
                        hier::IntVector::getZero(d_dim)));
                
                double* dudy = data_velocity_derivatives->getPointer(0);
                double* dvdx = data_velocity_derivatives->getPointer(1);
                
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
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                        new DerivativeFirstOrder(
                            "first order derivative in x-direction",
                            d_dim,
                            DIRECTION::X_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                        new DerivativeFirstOrder(
                            "first order derivative in y-direction",
                            d_dim,
                            DIRECTION::Y_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    // Compute dudy.
                    derivative_first_order_y->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[1],
                        patch_visible_box,
                        0,
                        0);
                    
                    // Compute dvdx.
                    derivative_first_order_x->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[0],
                        patch_visible_box,
                        1,
                        1);
                    
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
                            const int idx = (relative_idx_lo_0 + i) +
                                (relative_idx_lo_1 + j)*patch_dim_0;
                            
                            const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density;
                            
                            const double omega = dvdx[idx] - dudy[idx];
                            Omega_to_add += rho[idx_density]*omega*omega/((double) n_overlapped);
                        }
                    }
                }
                
                Omega_to_add = Omega_to_add*dx[0]*dx[1];
                Omega_integrated_local += Omega_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
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
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const HAMERS_SHARED_PTR<hier::Patch> patch = *ip;
                
                /*
                 * Get the patch geometry.
                 */
                
                const hier::Box& patch_box = patch->getBox();
                
                const hier::IntVector patch_dims = patch_box.numberCells();
                
                const int patch_dim_0 = patch_dims[0];
                const int patch_dim_1 = patch_dims[1];
                
                const hier::Index& patch_index_lo = patch_box.lower();
                
                const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                    HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch->getPatchGeometry()));
                
                const double* const dx = patch_geom->getDx();
                
                /*
                 * Register the patch, density and velocity in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", num_ghosts));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", num_ghosts));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->allocateMemoryForDerivedCellData();
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_density = flow_model_tmp->getCellData("DENSITY");
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity = flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                
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
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_ghosts_velocity >= hier::IntVector::getOne(d_dim)*d_num_ghosts_derivative);
#endif
                
                double Omega_to_add = 0.0;
                
                /*
                 * Initialize cell data for velocity derivatives and get pointers to the derivatives.
                 */
                
                HAMERS_SHARED_PTR<pdat::CellData<double> > data_velocity_derivatives(
                    new pdat::CellData<double>(patch_box, d_dim.getValue()*d_dim.getValue() - d_dim.getValue(),
                        hier::IntVector::getZero(d_dim)));
                
                double* dudy = data_velocity_derivatives->getPointer(0);
                double* dudz = data_velocity_derivatives->getPointer(1);
                double* dvdx = data_velocity_derivatives->getPointer(2);
                double* dvdz = data_velocity_derivatives->getPointer(3);
                double* dwdx = data_velocity_derivatives->getPointer(4);
                double* dwdy = data_velocity_derivatives->getPointer(5);
                
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
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_x(
                        new DerivativeFirstOrder(
                            "first order derivative in x-direction",
                            d_dim,
                            DIRECTION::X_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_y(
                        new DerivativeFirstOrder(
                            "first order derivative in y-direction",
                            d_dim,
                            DIRECTION::Y_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    HAMERS_SHARED_PTR<DerivativeFirstOrder> derivative_first_order_z(
                        new DerivativeFirstOrder(
                            "first order derivative in z-direction",
                            d_dim,
                            DIRECTION::Z_DIRECTION,
                            d_num_ghosts_derivative));
                    
                    // Compute dudy.
                    derivative_first_order_y->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[1],
                        patch_visible_box,
                        0,
                        0);
                    
                    // Compute dudz.
                    derivative_first_order_z->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[2],
                        patch_visible_box,
                        1,
                        0);
                    
                    // Compute dvdx.
                    derivative_first_order_x->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[0],
                        patch_visible_box,
                        2,
                        1);
                    
                    // Compute dvdz.
                    derivative_first_order_z->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[2],
                        patch_visible_box,
                        3,
                        1);
                    
                    // Compute dwdx.
                    derivative_first_order_x->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[0],
                        patch_visible_box,
                        4,
                        2);
                    
                    // Compute dwdy.
                    derivative_first_order_y->computeDerivative(
                        data_velocity_derivatives,
                        data_velocity,
                        dx[1],
                        patch_visible_box,
                        5,
                        2);
                    
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
                                
                                // Compute the linear indices.
                                const int idx = (relative_idx_lo_0 + i) +
                                    (relative_idx_lo_1 + j)*patch_dim_0 +
                                    (relative_idx_lo_2 + k)*patch_dim_0*
                                        patch_dim_1;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const double omega_x = dwdy[idx] - dvdz[idx];
                                const double omega_y = dudz[idx] - dwdx[idx];
                                const double omega_z = dvdx[idx] - dudy[idx];
                                
                                Omega_to_add += rho[idx_density]*(
                                    omega_x*omega_x + omega_y*omega_y + omega_z*omega_z)/((double) n_overlapped);
                            }
                        }
                    }
                }
                
                Omega_to_add = Omega_to_add*dx[0]*dx[1]*dx[2];
                Omega_integrated_local += Omega_to_add;
                
                /*
                 * Unregister the patch and data of all registered derived cell variables in the flow model.
                 */
                
                flow_model_tmp->unregisterPatch();
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
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output number of cells to a file.
 */
void
SBIStatisticsUtilities::outputNumberOfCells(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
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
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    if (d_dim == tbox::Dimension(1))
    {
        double num_cells_local = double(0);
        double num_cells_global = double(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
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
        
        mpi.Reduce(
            &num_cells_local,
            &num_cells_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the number of cells (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << num_cells_global;
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double num_cells_local = double(0);
        double num_cells_global = double(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
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
         * Reduction to get the global number of cells.
         */
        
        mpi.Reduce(
            &num_cells_local,
            &num_cells_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the number of cells (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << num_cells_global;
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double num_cells_local = double(0);
        double num_cells_global = double(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
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
        
        mpi.Reduce(
            &num_cells_local,
            &num_cells_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the number of cells (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << num_cells_global;
        }
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output weighted number of cells to a file.
 */
void
SBIStatisticsUtilities::outputWeightedNumberOfCells(
    const std::string& stat_dump_filename,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
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
    
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    
    if (d_dim == tbox::Dimension(1))
    {
        double weighted_num_cells_local = double(0);
        double weighted_num_cells_global = double(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the refinement ratio from the current level to the coarest level.
             */
            
            hier::IntVector ratioCurrentLevelToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioCurrentLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            if (li == 0)
            {
                ratioCurrentLevelToCoarestLevel = hier::IntVector::getOne(d_dim);
            }
            
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
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
                
                weighted_num_cells_local += double(interior_dims[0])*double(ratioCurrentLevelToCoarestLevel[0]);
            }
        }
        
        /*
         * Reduction to get the global number of cells.
         */
        
        mpi.Reduce(
            &weighted_num_cells_local,
            &weighted_num_cells_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the number of cells (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << weighted_num_cells_global;
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        double weighted_num_cells_local = double(0);
        double weighted_num_cells_global = double(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the refinement ratio from the current level to the coarest level.
             */
            
            hier::IntVector ratioCurrentLevelToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioCurrentLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            if (li == 0)
            {
                ratioCurrentLevelToCoarestLevel = hier::IntVector::getOne(d_dim);
            }
            
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
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
                    double(ratioCurrentLevelToCoarestLevel[0]);
            }
        }
        
        /*
         * Reduction to get the global number of cells.
         */
        
        mpi.Reduce(
            &weighted_num_cells_local,
            &weighted_num_cells_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the number of cells (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << weighted_num_cells_global;
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        double weighted_num_cells_local = double(0);
        double weighted_num_cells_global = double(0);
        
        for (int li = 0; li < num_levels; li++)
        {
            /*
             * Get the refinement ratio from the current level to the coarest level.
             */
            
            hier::IntVector ratioCurrentLevelToCoarestLevel =
                patch_hierarchy->getRatioToCoarserLevel(li);
            for (int lii = li - 1; lii > 0 ; lii--)
            {
                ratioCurrentLevelToCoarestLevel *= patch_hierarchy->getRatioToCoarserLevel(lii);
            }
            
            if (li == 0)
            {
                ratioCurrentLevelToCoarestLevel = hier::IntVector::getOne(d_dim);
            }
            
            /*
             * Get the current patch level.
             */
            
            HAMERS_SHARED_PTR<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
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
                
                weighted_num_cells_local += double(interior_dims[0])*double(interior_dims[1])*double(interior_dims[2])*
                    double(ratioCurrentLevelToCoarestLevel[0]);
            }
        }
        
        /*
         * Reduction to get the global number of cells.
         */
        
        mpi.Reduce(
            &weighted_num_cells_local,
            &weighted_num_cells_global,
            1,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Output the number of cells (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << weighted_num_cells_global;
        }
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Compute statisitcal quantities.
 */
void
FlowModelStatisticsUtilitiesFourEqnConservative::computeStatisticalQuantities(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double statistics_data_time)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(data_context);
    NULL_USE(statistics_data_time);
}


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
            else if (statistical_quantity_key == "CENTROID_X")
            {
                f_out << "\t" << "CENTROID_X           ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_X")
            {
                f_out << "\t" << "INTERFACE_MIN_X      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_X")
            {
                f_out << "\t" << "INTERFACE_MAX_X      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_Y")
            {
                f_out << "\t" << "INTERFACE_MIN_Y      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_Y")
            {
                f_out << "\t" << "INTERFACE_MAX_Y      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MIN_Z")
            {
                f_out << "\t" << "INTERFACE_MIN_Z      ";
            }
            else if (statistical_quantity_key == "INTERFACE_MAX_Z")
            {
                f_out << "\t" << "INTERFACE_MAX_Z      ";
            }
            else if (statistical_quantity_key == "CIRCULATION")
            {
                f_out << "\t" << "CIRCULATION          ";
            }
            else if (statistical_quantity_key == "ENSTROPHY_INT")
            {
                f_out << "\t" << "ENSTROPHY_INT        ";
            }
            else if (statistical_quantity_key == "SCAL_DISS_RAT_INT")
            {
                f_out << "\t" << "SCAL_DISS_RAT_INT    ";
            }
            else if (statistical_quantity_key == "NUM_CELLS")
            {
                f_out << "\t" << "NUM_CELLS            ";
            }
            else if (statistical_quantity_key == "WEIGHTED_NUM_CELLS")
            {
                f_out << "\t" << "WEIGHTED_NUM_CELLS   ";
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
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
    const double output_time)
{
    NULL_USE(output_time);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    HAMERS_SHARED_PTR<SBIStatisticsUtilities> sbi_statistics_utilities(
        new SBIStatisticsUtilities(
            "SBI statistics utilities",
            d_dim,
            d_grid_geometry,
            d_num_species,
            d_flow_model,
            d_equation_of_state_mixing_rules,
            d_equation_of_mass_diffusivity_mixing_rules,
            d_equation_of_shear_viscosity_mixing_rules,
            d_equation_of_bulk_viscosity_mixing_rules,
            d_equation_of_thermal_conductivity_mixing_rules));
    
    // Loop over statistical quantities.
    for (int qi = 0; qi < static_cast<int>(d_statistical_quantities.size()); qi++)
    {
        // Get the key of the current variable.
        std::string statistical_quantity_key = d_statistical_quantities[qi];
        
        if (statistical_quantity_key == "MIXING_WIDTH_X")
        {
            sbi_statistics_utilities->outputMixingWidthInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "CENTROID_X")
        {
            sbi_statistics_utilities->outputCentroidInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_X")
        {
            sbi_statistics_utilities->outputInterfaceMinInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_X")
        {
            sbi_statistics_utilities->outputInterfaceMaxInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_Y")
        {
            sbi_statistics_utilities->outputInterfaceMinInYDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_Y")
        {
            sbi_statistics_utilities->outputInterfaceMaxInYDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MIN_Z")
        {
            sbi_statistics_utilities->outputInterfaceMinInZDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "INTERFACE_MAX_Z")
        {
            sbi_statistics_utilities->outputInterfaceMaxInZDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "CIRCULATION")
        {
            sbi_statistics_utilities->outputCirculation(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "ENSTROPHY_INT")
        {
            sbi_statistics_utilities->outputEnstrophyIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "SCAL_DISS_RAT_INT")
        {
            sbi_statistics_utilities->outputScalarDissipationRateIntegrated(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "NUM_CELLS")
        {
            sbi_statistics_utilities->outputNumberOfCells(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "WEIGHTED_NUM_CELLS")
        {
            sbi_statistics_utilities->outputWeightedNumberOfCells(
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
                << "' found."
                << std::endl);
        }
    }
}
