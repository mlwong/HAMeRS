#include "flow/flow_models/four-eqn_conservative/FlowModelStatisticsUtilitiesFourEqnConservative.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"

#include <fstream>

class SecondMomentStatisticsUtilities
{
    public:
        SecondMomentStatisticsUtilities(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::weak_ptr<FlowModel> flow_model,
            const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules,
            const boost::shared_ptr<EquationOfMassDiffusivityMixingRules> equation_of_mass_diffusivity_mixing_rules,
            const boost::shared_ptr<EquationOfShearViscosityMixingRules> equation_of_shear_viscosity_mixing_rules,
            const boost::shared_ptr<EquationOfBulkViscosityMixingRules> equation_of_bulk_viscosity_mixing_rules,
            const boost::shared_ptr<EquationOfThermalConductivityMixingRules> equation_of_thermal_conductivity_mixing_rules):
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
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output mixedness in x-direction to a file.
         */
        void
        outputMixednessInXDirection(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output integrated or mean velocity of turbulent mass flux in x-direction inside mixing layer with
         * assumed homogeneity in yz-plane to a file.
         */
        void
        outputTurbulentMassFluxVelocityInXDirectionInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean turbulent mass flux in x-direction inside mixing layer with assumed
         * homogeneity in yz-plane to a file.
         */
        void
        outputTurbulentMassFluxInXDirectionInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean density specific volume covariance inside mixing layer with assumed
         * homogeneity in yz-plane to a file.
         */
        void
        outputDensitySpecificVolumeCovarianceInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean density multiplied by density specific volume covariance inside mixing layer
         * with assumed homogeneity in yz-plane to a file.
         */
        void
        outputDensityTimesDensitySpecificVolumeCovarianceInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean Reynolds normal stress component in x-direction inside mixing layer with
         * assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsNormalStressInXDirectionInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean density multiplied by Reynolds normal stress component in x-direction
         * inside mixing layer with assumed homogeneity in yz-plane to a file.
         */
        void
        outputDensityTimesReynoldsNormalStressInXDirectionInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean Reynolds normal stress component in y-direction inside mixing layer with
         * assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsNormalStressInYDirectionInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean density multiplied by Reynolds normal stress component in y-direction
         * inside mixing layer with assumed homogeneity in yz-plane to a file.
         */
        void
        outputDensityTimesReynoldsNormalStressInYDirectionInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean Reynolds normal stress component in z-direction inside mixing layer with
         * assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsNormalStressInZDirectionInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean density multiplied by Reynolds normal stress component in z-direction
         * inside mixing layer with assumed homogeneity in yz-plane to a file.
         */
        void
        outputDensityTimesReynoldsNormalStressInZDirectionInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean Reynolds shear stress component in x- and y-directions inside mixing layer
         * with assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsShearStressInXYDirectionsInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean density multiplied by Reynolds shear stress component in x- and y-directions
         * inside mixing layer with assumed homogeneity in yz-plane to a file.
         */
        void
        outputDensityTimesReynoldsShearStressInXYDirectionsInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean Reynolds shear stress component in x- and z-directions inside mixing layer
         * with assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsShearStressInXZDirectionsInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean density multiplied by Reynolds shear stress component in x- and z-directions
         * inside mixing layer with assumed homogeneity in yz-plane to a file.
         */
        void
        outputDensityTimesReynoldsShearStressInXZDirectionsInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean Reynolds shear stress component in y- and z-directions inside mixing layer
         * with assumed homogeneity in yz-plane to a file.
         */
        void
        outputReynoldsShearStressInYZDirectionsInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output integrated or mean density multiplied by Reynolds shear stress component in y- and z-directions
         * inside mixing layer with assumed homogeneity in yz-plane to a file.
         */
        void
        outputDensityTimesReynoldsShearStressInYZDirectionsInMixingLayerWithHomogeneityInYZPlane(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool integrated);
        
        /*
         * Output number of cells to a file.
         */
        void
        outputNumberOfCells(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Output weighted number of cells to a file.
         */
        void
        outputWeightedNumberOfCells(
            const std::string& stat_dump_filename,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
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
         * boost::shared_ptr to the grid geometry.
         */
        const boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * boost::weak_ptr to FlowModel.
         */
        const boost::weak_ptr<FlowModel> d_flow_model;
        
        /*
         * boost::shared_ptr to EquationOfStateMixingRules.
         */
        const boost::shared_ptr<EquationOfStateMixingRules>
            d_equation_of_state_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfMassDiffusivityMixingRules.
         */
        const boost::shared_ptr<EquationOfMassDiffusivityMixingRules>
            d_equation_of_mass_diffusivity_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfShearViscosityMixingRules.
         */
        const boost::shared_ptr<EquationOfShearViscosityMixingRules>
            d_equation_of_shear_viscosity_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfBulkViscosityMixingRules.
         */
        const boost::shared_ptr<EquationOfBulkViscosityMixingRules>
            d_equation_of_bulk_viscosity_mixing_rules;
        
        /*
         * boost::shared_ptr to EquationOfThermalConductivityMixingRules.
         */
        const boost::shared_ptr<EquationOfThermalConductivityMixingRules>
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
SecondMomentStatisticsUtilities::outputMixingWidthInXDirection(
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
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
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
                 * Register the patch and mole fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mole fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                double* X = data_mole_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                
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
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_mole_fractions;
                        
                        const double value_to_add = fmax(fmin(X[idx], 1.0), 0.0)/((double) n_overlapped);
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            X_avg_local[idx_fine] += value_to_add;
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
            X_avg_local,
            X_avg_global,
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
                W += X_avg_global[i]*(1.0 - X_avg_global[i]);
            }
            
            const double dx_finest = L_x/finest_level_dim_0;
            
            W = 4.0*W*dx_finest;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << W;
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_x = x_hi[0] - x_lo[0];
        const double L_y = x_hi[1] - x_lo[1];
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
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
                 * Register the patch and mole fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mole fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                double* X = data_mole_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                
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
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions;
                            
                            const double value_to_add = fmax(fmin(X[idx], 1.0), 0.0)*weight/((double) n_overlapped);
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                X_avg_local[idx_fine] += value_to_add;
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
            X_avg_local,
            X_avg_global,
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
                W += X_avg_global[i]*(1.0 - X_avg_global[i]);
            }
            
            const double dx_finest = L_x/finest_level_dim_0;
            
            W = 4.0*W*dx_finest;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << W;
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
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
                 * Register the patch and mole fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mole fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                double* X = data_mole_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const double value_to_add = fmax(fmin(X[idx], 1.0), 0.0)*weight/((double) n_overlapped);
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += value_to_add;
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
            X_avg_local,
            X_avg_global,
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
                W += X_avg_global[i]*(1.0 - X_avg_global[i]);
            }
            
            const double dx_finest = L_x/finest_level_dim_0;
            
            W = 4.0*W*dx_finest;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << W;
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output mixedness in x-direction to a file.
 */
void
SecondMomentStatisticsUtilities::outputMixednessInXDirection(
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
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_product_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_product_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            X_product_avg_local[i] = 0.0;
            X_product_avg_global[i] = 0.0;
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
                 * Register the patch and mole fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mole fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                double* X = data_mole_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                
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
                        
                        const int idx = relative_idx_lo_0 + i + num_ghosts_0_mole_fractions;
                        
                        const double X_bounded = fmax(fmin(X[idx], 1.0), 0.0);
                        
                        const double weight_local = 1.0/((double) n_overlapped);
                        
                        const double value_to_add = X_bounded*weight_local;
                        const double product_to_add = X_bounded*(1.0 - X_bounded)*weight_local;
                        
                        for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                        {
                            const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                            
                            X_avg_local[idx_fine] += value_to_add;
                            X_product_avg_local[idx_fine] += product_to_add;
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
         * Reduction to get the global averages.
         */
        
        mpi.Reduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            X_product_avg_local,
            X_product_avg_global,
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
                num += X_product_avg_global[i];
            }
            
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                den += X_avg_global[i]*(1.0 - X_avg_global[i]);
            }
            
            const double Theta = num/den;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Theta;
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(X_product_avg_local);
        std::free(X_product_avg_global);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_product_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_product_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            X_product_avg_local[i] = 0.0;
            X_product_avg_global[i] = 0.0;
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
                 * Register the patch and mole fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mole fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                double* X = data_mole_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                
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
                            
                            const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions;
                            
                            const double X_bounded = fmax(fmin(X[idx], 1.0), 0.0);
                            
                            const double weight_local = weight/((double) n_overlapped);
                            
                            const double value_to_add = X_bounded*weight_local;
                            const double product_to_add = X_bounded*(1.0 - X_bounded)*weight_local;
                            
                            for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                            {
                                const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                
                                X_avg_local[idx_fine] += value_to_add;
                                X_product_avg_local[idx_fine] += product_to_add;
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
         * Reduction to get the global averages.
         */
        
        mpi.Reduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            X_product_avg_local,
            X_product_avg_global,
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
                num += X_product_avg_global[i];
            }
            
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                den += X_avg_global[i]*(1.0 - X_avg_global[i]);
            }
            
            const double Theta = num/den;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Theta;
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(X_product_avg_local);
        std::free(X_product_avg_global);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int finest_level_dim_0 = finest_level_dims[0];
        
        /*
         * Get the size of the physical domain.
         */
        
        const double L_y = x_hi[1] - x_lo[1];
        const double L_z = x_hi[2] - x_lo[2];
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_product_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_product_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            X_product_avg_local[i] = 0.0;
            X_product_avg_global[i] = 0.0;
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
                 * Register the patch and mole fractions in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointer to first mole fraction data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                double* X = data_mole_fractions->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const double X_bounded = fmax(fmin(X[idx], 1.0), 0.0);
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double value_to_add = X_bounded*weight_local;
                                const double product_to_add = X_bounded*(1.0 - X_bounded)*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += value_to_add;
                                    X_product_avg_local[idx_fine] += product_to_add;
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
         * Reduction to get the global averages.
         */
        
        mpi.Reduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        mpi.Reduce(
            X_product_avg_local,
            X_product_avg_global,
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
                num += X_product_avg_global[i];
            }
            
            for (int i = 0; i < finest_level_dim_0; i++)
            {
                den += X_avg_global[i]*(1.0 - X_avg_global[i]);
            }
            
            const double Theta = num/den;
            
            f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                  << "\t" << Theta;
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(X_product_avg_local);
        std::free(X_product_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean velocity of turbulent mass flux in x-direction inside mixing layer with
 * assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputTurbulentMassFluxVelocityInXDirectionInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'a1_INT_HOMO_YZ'/'a1_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'a1_INT_HOMO_YZ'/'a1_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'a1_INT_HOMO_YZ'/'a1_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* u_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
            u_avg_local[i] = 0.0;
            u_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double u_to_add = u[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    u_avg_local[idx_fine] += u_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        double* rho_p_u_p_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_p_u_p_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_p_u_p_avg_local[i] = 0.0;
            rho_p_u_p_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double rho_p = rho[idx_density] - rho_avg_global[idx_fine];
                                    const double u_p = u[idx_velocity] - u_avg_global[idx_fine];
                                    
                                    rho_p_u_p_avg_local[idx_fine] += rho_p*u_p*weight_local;
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
        
        mpi.Reduce(
            rho_p_u_p_avg_local,
            rho_p_u_p_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean velocity of turbulent mass flux inside mixing layer
         * (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double a_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    a_sum += rho_p_u_p_avg_global[i]*dx/rho_avg_global[i];
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << a_sum;
            }
            else
            {
                double a_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        a_sum += rho_p_u_p_avg_global[i]/rho_avg_global[i];
                        count++;
                    }
                }
                
                const double a_mean = a_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << a_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(u_avg_local);
        std::free(u_avg_global);
        std::free(rho_p_u_p_avg_local);
        std::free(rho_p_u_p_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean turbulent mass flux in x-direction inside mixing layer with assumed
 * homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputTurbulentMassFluxInXDirectionInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'ra1_INT_HOMO_YZ'/'ra1_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'ra1_INT_HOMO_YZ'/'ra1_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'ra1_INT_HOMO_YZ'/'ra1_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* u_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
            u_avg_local[i] = 0.0;
            u_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double u_to_add = u[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    u_avg_local[idx_fine] += u_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            u_avg_local,
            u_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        double* rho_p_u_p_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_p_u_p_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_p_u_p_avg_local[i] = 0.0;
            rho_p_u_p_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double rho_p = rho[idx_density] - rho_avg_global[idx_fine];
                                    const double u_p = u[idx_velocity] - u_avg_global[idx_fine];
                                    
                                    rho_p_u_p_avg_local[idx_fine] += rho_p*u_p*weight_local;
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
        
        mpi.Reduce(
            rho_p_u_p_avg_local,
            rho_p_u_p_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean turbulent mass flux inside mixing layer
         * (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double rho_a_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    rho_a_sum += rho_p_u_p_avg_global[i]*dx;
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_a_sum;
            }
            else
            {
                double rho_a_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        rho_a_sum += rho_p_u_p_avg_global[i];
                        count++;
                    }
                }
                
                const double rho_a_mean = rho_a_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_a_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(u_avg_local);
        std::free(u_avg_global);
        std::free(rho_p_u_p_avg_local);
        std::free(rho_p_u_p_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean density specific volume covariance inside mixing layer with assumed
 * homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputDensitySpecificVolumeCovarianceInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'b_INT_HOMO_YZ'/'b_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'b_INT_HOMO_YZ'/'b_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'b_INT_HOMO_YZ'/'b_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* v_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* v_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
            v_avg_local[i] = 0.0;
            v_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions and density in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions and density data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double v_to_add = 1.0/rho[idx_density]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    v_avg_local[idx_fine] += v_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            v_avg_local,
            v_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        double* rho_p_v_p_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_p_v_p_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_p_v_p_avg_local[i] = 0.0;
            rho_p_v_p_avg_global[i] = 0.0;
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
                 * Register the patch and density in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
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
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
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
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double rho_p = rho[idx_density] - rho_avg_global[idx_fine];
                                    const double v_p = 1.0/rho[idx_density] - v_avg_global[idx_fine];
                                    
                                    rho_p_v_p_avg_local[idx_fine] += rho_p*v_p*weight_local;
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
        
        mpi.Reduce(
            rho_p_v_p_avg_local,
            rho_p_v_p_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean density specific volume covariance inside mixing layer
         * (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double b_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    b_sum += -rho_p_v_p_avg_global[i]*dx;
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << b_sum;
            }
            else
            {
                double b_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        b_sum += -rho_p_v_p_avg_global[i];
                        count++;
                    }
                }
                
                const double b_mean = b_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << b_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(v_avg_local);
        std::free(v_avg_global);
        std::free(rho_p_v_p_avg_local);
        std::free(rho_p_v_p_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean density multiplied by density specific volume covariance inside mixing layer
 * with assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputDensityTimesDensitySpecificVolumeCovarianceInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'rb_INT_HOMO_YZ'/'rb_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'rb_INT_HOMO_YZ'/'rb_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'rb_INT_HOMO_YZ'/'rb_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* v_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* v_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
            v_avg_local[i] = 0.0;
            v_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions and density in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions and density data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double v_to_add = 1.0/rho[idx_density]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    v_avg_local[idx_fine] += v_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            v_avg_local,
            v_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        double* rho_p_v_p_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_p_v_p_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_p_v_p_avg_local[i] = 0.0;
            rho_p_v_p_avg_global[i] = 0.0;
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
                 * Register the patch and density in the flow model and compute the
                 * corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
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
                
                const int num_ghosts_0_density = num_ghosts_density[0];
                const int num_ghosts_1_density = num_ghosts_density[1];
                const int num_ghosts_2_density = num_ghosts_density[2];
                const int ghostcell_dim_0_density = ghostcell_dims_density[0];
                const int ghostcell_dim_1_density = ghostcell_dims_density[1];
                
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
                                 * Compute the linear index and the data to add.
                                 */
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double rho_p = rho[idx_density] - rho_avg_global[idx_fine];
                                    const double v_p = 1.0/rho[idx_density] - v_avg_global[idx_fine];
                                    
                                    rho_p_v_p_avg_local[idx_fine] += rho_p*v_p*weight_local;
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
        
        mpi.Reduce(
            rho_p_v_p_avg_local,
            rho_p_v_p_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean density multiplied by density specific volume covariance
         * inside mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double rb_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    rb_sum += -rho_avg_global[i]*rho_p_v_p_avg_global[i]*dx;
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rb_sum;
            }
            else
            {
                double rb_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        rb_sum += -rho_avg_global[i]*rho_p_v_p_avg_global[i];
                        count++;
                    }
                }
                
                const double rb_mean = rb_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rb_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(v_avg_local);
        std::free(v_avg_global);
        std::free(rho_p_v_p_avg_local);
        std::free(rho_p_v_p_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean Reynolds normal stress component in x-direction inside mixing layer with
 * assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputReynoldsNormalStressInXDirectionInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'R11_INT_HOMO_YZ'/'R11_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'R11_INT_HOMO_YZ'/'R11_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'R11_INT_HOMO_YZ'/'R11_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
            rho_u_avg_local[i] = 0.0;
            rho_u_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_u_to_add = rho[idx_density]*u[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_u_avg_local[idx_fine] += rho_u_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
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
        
        double* rho_u_pp_u_pp_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_pp_u_pp_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_u_pp_u_pp_avg_local[i] = 0.0;
            rho_u_pp_u_pp_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double u_pp = u[idx_velocity] - rho_u_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    rho_u_pp_u_pp_avg_local[idx_fine] +=
                                        rho[idx_density]*u_pp*u_pp*weight_local;
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
        
        mpi.Reduce(
            rho_u_pp_u_pp_avg_local,
            rho_u_pp_u_pp_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean Reynolds normal stress component in x-direction inside
         * mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double R11_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    R11_sum += rho_u_pp_u_pp_avg_global[i]*dx/rho_avg_global[i];
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << R11_sum;
            }
            else
            {
                double R11_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        R11_sum += rho_u_pp_u_pp_avg_global[i]/rho_avg_global[i];
                        count++;
                    }
                }
                
                const double R11_mean = R11_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << R11_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_u_avg_local);
        std::free(rho_u_avg_global);
        std::free(rho_u_pp_u_pp_avg_local);
        std::free(rho_u_pp_u_pp_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean density multiplied by Reynolds normal stress component in x-direction
 * inside mixing layer with assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputDensityTimesReynoldsNormalStressInXDirectionInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'rR11_INT_HOMO_YZ'/'rR11_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'rR11_INT_HOMO_YZ'/'rR11_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'rR11_INT_HOMO_YZ'/'rR11_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
            rho_u_avg_local[i] = 0.0;
            rho_u_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_u_to_add = rho[idx_density]*u[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_u_avg_local[idx_fine] += rho_u_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
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
        
        double* rho_u_pp_u_pp_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_pp_u_pp_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_u_pp_u_pp_avg_local[i] = 0.0;
            rho_u_pp_u_pp_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double u_pp = u[idx_velocity] - rho_u_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    rho_u_pp_u_pp_avg_local[idx_fine] +=
                                        rho[idx_density]*u_pp*u_pp*weight_local;
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
        
        mpi.Reduce(
            rho_u_pp_u_pp_avg_local,
            rho_u_pp_u_pp_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean density multiplied by Reynolds normal stress component
         * in x-direction inside mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double rho_R11_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    rho_R11_sum += rho_u_pp_u_pp_avg_global[i]*dx;
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_R11_sum;
            }
            else
            {
                double rho_R11_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        rho_R11_sum += rho_u_pp_u_pp_avg_global[i];
                        count++;
                    }
                }
                
                const double rho_R11_mean = rho_R11_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_R11_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_u_avg_local);
        std::free(rho_u_avg_global);
        std::free(rho_u_pp_u_pp_avg_local);
        std::free(rho_u_pp_u_pp_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean Reynolds normal stress component in y-direction inside mixing layer with
 * assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputReynoldsNormalStressInYDirectionInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'R22_INT_HOMO_YZ'/'R22_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'R22_INT_HOMO_YZ'/'R22_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'R22_INT_HOMO_YZ'/'R22_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
                double* v = data_velocity->getPointer(1);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_v_to_add = rho[idx_density]*v[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_v_avg_local[idx_fine] += rho_v_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_v_avg_local,
            rho_v_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        double* rho_v_pp_v_pp_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_pp_v_pp_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_v_pp_v_pp_avg_local[i] = 0.0;
            rho_v_pp_v_pp_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double v_pp = v[idx_velocity] - rho_v_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    rho_v_pp_v_pp_avg_local[idx_fine] +=
                                        rho[idx_density]*v_pp*v_pp*weight_local;
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
        
        mpi.Reduce(
            rho_v_pp_v_pp_avg_local,
            rho_v_pp_v_pp_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean Reynolds normal stress component in y-direction inside
         * mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double R22_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    R22_sum += rho_v_pp_v_pp_avg_global[i]*dx/rho_avg_global[i];
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << R22_sum;
            }
            else
            {
                double R22_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        R22_sum += rho_v_pp_v_pp_avg_global[i]/rho_avg_global[i];
                        count++;
                    }
                }
                
                const double R22_mean = R22_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << R22_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_v_avg_local);
        std::free(rho_v_avg_global);
        std::free(rho_v_pp_v_pp_avg_local);
        std::free(rho_v_pp_v_pp_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean density multiplied by Reynolds normal stress component in y-direction
 * inside mixing layer with assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputDensityTimesReynoldsNormalStressInYDirectionInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'rR22_INT_HOMO_YZ'/'rR22_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'rR22_INT_HOMO_YZ'/'rR22_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'rR22_INT_HOMO_YZ'/'rR22_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
                double* v = data_velocity->getPointer(1);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_v_to_add = rho[idx_density]*v[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_v_avg_local[idx_fine] += rho_v_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_v_avg_local,
            rho_v_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        double* rho_v_pp_v_pp_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_pp_v_pp_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_v_pp_v_pp_avg_local[i] = 0.0;
            rho_v_pp_v_pp_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double v_pp = v[idx_velocity] - rho_v_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    rho_v_pp_v_pp_avg_local[idx_fine] +=
                                        rho[idx_density]*v_pp*v_pp*weight_local;
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
        
        mpi.Reduce(
            rho_v_pp_v_pp_avg_local,
            rho_v_pp_v_pp_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean density multiplied by Reynolds normal stress component
         * in y-direction inside mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double rho_R22_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    rho_R22_sum += rho_v_pp_v_pp_avg_global[i]*dx;
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_R22_sum;
            }
            else
            {
                double rho_R22_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        rho_R22_sum += rho_v_pp_v_pp_avg_global[i];
                        count++;
                    }
                }
                
                const double rho_R22_mean = rho_R22_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_R22_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_v_avg_local);
        std::free(rho_v_avg_global);
        std::free(rho_v_pp_v_pp_avg_local);
        std::free(rho_v_pp_v_pp_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean Reynolds normal stress component in z-direction inside mixing layer with
 * assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputReynoldsNormalStressInZDirectionInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'R33_INT_HOMO_YZ'/'R33_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'R33_INT_HOMO_YZ'/'R33_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'R33_INT_HOMO_YZ'/'R33_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
                double* w = data_velocity->getPointer(2);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_w_to_add = rho[idx_density]*w[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_w_avg_local[idx_fine] += rho_w_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_w_avg_local,
            rho_w_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        double* rho_w_pp_w_pp_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_pp_w_pp_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_w_pp_w_pp_avg_local[i] = 0.0;
            rho_w_pp_w_pp_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double w_pp = w[idx_velocity] - rho_w_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    rho_w_pp_w_pp_avg_local[idx_fine] +=
                                        rho[idx_density]*w_pp*w_pp*weight_local;
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
        
        mpi.Reduce(
            rho_w_pp_w_pp_avg_local,
            rho_w_pp_w_pp_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean Reynolds normal stress component in z-direction inside
         * mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double R33_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    R33_sum += rho_w_pp_w_pp_avg_global[i]*dx/rho_avg_global[i];
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << R33_sum;
            }
            else
            {
                double R33_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        R33_sum += rho_w_pp_w_pp_avg_global[i]/rho_avg_global[i];
                        count++;
                    }
                }
                
                const double R33_mean = R33_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << R33_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_w_avg_local);
        std::free(rho_w_avg_global);
        std::free(rho_w_pp_w_pp_avg_local);
        std::free(rho_w_pp_w_pp_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean density multiplied by Reynolds normal stress component in z-direction
 * inside mixing layer with assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputDensityTimesReynoldsNormalStressInZDirectionInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'rR33_INT_HOMO_YZ'/'rR33_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'rR33_INT_HOMO_YZ'/'rR33_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'rR33_INT_HOMO_YZ'/'rR33_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
                double* w = data_velocity->getPointer(2);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_w_to_add = rho[idx_density]*w[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_w_avg_local[idx_fine] += rho_w_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_w_avg_local,
            rho_w_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        double* rho_w_pp_w_pp_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_pp_w_pp_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_w_pp_w_pp_avg_local[i] = 0.0;
            rho_w_pp_w_pp_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double w_pp = w[idx_velocity] - rho_w_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    rho_w_pp_w_pp_avg_local[idx_fine] +=
                                        rho[idx_density]*w_pp*w_pp*weight_local;
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
        
        mpi.Reduce(
            rho_w_pp_w_pp_avg_local,
            rho_w_pp_w_pp_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean density multiplied Reynolds normal stress component
         * in z-direction inside mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double rho_R33_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    rho_R33_sum += rho_w_pp_w_pp_avg_global[i]*dx;
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_R33_sum;
            }
            else
            {
                double rho_R33_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        rho_R33_sum += rho_w_pp_w_pp_avg_global[i];
                        count++;
                    }
                }
                
                const double rho_R33_mean = rho_R33_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_R33_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_w_avg_local);
        std::free(rho_w_avg_global);
        std::free(rho_w_pp_w_pp_avg_local);
        std::free(rho_w_pp_w_pp_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean Reynolds shear stress component in x- and y-directions inside mixing layer
 * with assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputReynoldsShearStressInXYDirectionsInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'R12_INT_HOMO_YZ'/'R12_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'R12_INT_HOMO_YZ'/'R12_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'R12_INT_HOMO_YZ'/'R12_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
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
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_u_to_add = rho[idx_density]*u[idx_velocity]*weight_local;
                                const double rho_v_to_add = rho[idx_density]*v[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_u_avg_local[idx_fine] += rho_u_to_add;
                                    rho_v_avg_local[idx_fine] += rho_v_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
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
        
        double* rho_u_pp_v_pp_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_pp_v_pp_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_u_pp_v_pp_avg_local[i] = 0.0;
            rho_u_pp_v_pp_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double u_pp = u[idx_velocity] - rho_u_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    const double v_pp = v[idx_velocity] - rho_v_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    rho_u_pp_v_pp_avg_local[idx_fine] +=
                                        rho[idx_density]*u_pp*v_pp*weight_local;
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
        
        mpi.Reduce(
            rho_u_pp_v_pp_avg_local,
            rho_u_pp_v_pp_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean Reynolds shear stress component in x- and y-directions
         * inside mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double R12_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                        R12_sum += rho_u_pp_v_pp_avg_global[i]*dx/rho_avg_global[i];
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << R12_sum;
            }
            else
            {
                double R12_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        R12_sum += rho_u_pp_v_pp_avg_global[i]/rho_avg_global[i];
                        count++;
                    }
                }
                
                const double R12_mean = R12_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << R12_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_u_avg_local);
        std::free(rho_u_avg_global);
        std::free(rho_v_avg_local);
        std::free(rho_v_avg_global);
        std::free(rho_u_pp_v_pp_avg_local);
        std::free(rho_u_pp_v_pp_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean density multiplied by Reynolds shear stress component in x- and y-directions
 * inside mixing layer with assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputDensityTimesReynoldsShearStressInXYDirectionsInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'rR12_INT_HOMO_YZ'/'rR12_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'rR12_INT_HOMO_YZ'/'rR12_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'rR12_INT_HOMO_YZ'/'rR12_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
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
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_u_to_add = rho[idx_density]*u[idx_velocity]*weight_local;
                                const double rho_v_to_add = rho[idx_density]*v[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_u_avg_local[idx_fine] += rho_u_to_add;
                                    rho_v_avg_local[idx_fine] += rho_v_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
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
        
        double* rho_u_pp_v_pp_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_pp_v_pp_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_u_pp_v_pp_avg_local[i] = 0.0;
            rho_u_pp_v_pp_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double u_pp = u[idx_velocity] - rho_u_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    const double v_pp = v[idx_velocity] - rho_v_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    rho_u_pp_v_pp_avg_local[idx_fine] +=
                                        rho[idx_density]*u_pp*v_pp*weight_local;
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
        
        mpi.Reduce(
            rho_u_pp_v_pp_avg_local,
            rho_u_pp_v_pp_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean density multiplied by Reynolds shear stress component
         * in x- and y-directions inside mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double rho_R12_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                        rho_R12_sum += rho_u_pp_v_pp_avg_global[i]*dx;
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_R12_sum;
            }
            else
            {
                double rho_R12_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        rho_R12_sum += rho_u_pp_v_pp_avg_global[i];
                        count++;
                    }
                }
                
                const double rho_R12_mean = rho_R12_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_R12_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_u_avg_local);
        std::free(rho_u_avg_global);
        std::free(rho_v_avg_local);
        std::free(rho_v_avg_global);
        std::free(rho_u_pp_v_pp_avg_local);
        std::free(rho_u_pp_v_pp_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean Reynolds shear stress component in x- and z-directions inside mixing layer
 * with assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputReynoldsShearStressInXZDirectionsInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'R13_INT_HOMO_YZ'/'R13_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'R13_INT_HOMO_YZ'/'R13_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'R13_INT_HOMO_YZ'/'R13_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
            rho_u_avg_local[i] = 0.0;
            rho_u_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                double* w = data_velocity->getPointer(2);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_u_to_add = rho[idx_density]*u[idx_velocity]*weight_local;
                                const double rho_w_to_add = rho[idx_density]*w[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_u_avg_local[idx_fine] += rho_u_to_add;
                                    rho_w_avg_local[idx_fine] += rho_w_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
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
            rho_w_avg_local,
            rho_w_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        double* rho_u_pp_w_pp_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_pp_w_pp_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_u_pp_w_pp_avg_local[i] = 0.0;
            rho_u_pp_w_pp_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double u_pp = u[idx_velocity] - rho_u_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    const double w_pp = w[idx_velocity] - rho_w_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    rho_u_pp_w_pp_avg_local[idx_fine] +=
                                        rho[idx_density]*u_pp*w_pp*weight_local;
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
        
        mpi.Reduce(
            rho_u_pp_w_pp_avg_local,
            rho_u_pp_w_pp_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean Reynolds shear stress component in x- and z-directions
         * inside mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double R13_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    R13_sum += rho_u_pp_w_pp_avg_global[i]*dx/rho_avg_global[i];
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << R13_sum;
            }
            else
            {
                double R13_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        R13_sum += rho_u_pp_w_pp_avg_global[i]/rho_avg_global[i];
                        count++;
                    }
                }
                
                const double R13_mean = R13_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << R13_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_u_avg_local);
        std::free(rho_u_avg_global);
        std::free(rho_w_avg_local);
        std::free(rho_w_avg_global);
        std::free(rho_u_pp_w_pp_avg_local);
        std::free(rho_u_pp_w_pp_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean density multiplied by Reynolds shear stress component in x- and z-directions
 * inside mixing layer with assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputDensityTimesReynoldsShearStressInXZDirectionsInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'rR13_INT_HOMO_YZ'/'rR13_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'rR13_INT_HOMO_YZ'/'rR13_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'rR13_INT_HOMO_YZ'/'rR13_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
            rho_u_avg_local[i] = 0.0;
            rho_u_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
                double* w = data_velocity->getPointer(2);
                
                const hier::BoxContainer& patch_visible_boxes =
                    flattened_hierarchy->getVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::BoxContainer& patch_overlapped_visible_boxes =
                    flattened_hierarchy->getOverlappedVisibleBoxes(
                        patch_box,
                        li);
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_u_to_add = rho[idx_density]*u[idx_velocity]*weight_local;
                                const double rho_w_to_add = rho[idx_density]*w[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
                                    rho_u_avg_local[idx_fine] += rho_u_to_add;
                                    rho_w_avg_local[idx_fine] += rho_w_to_add;
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
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
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
            rho_w_avg_local,
            rho_w_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        double* rho_u_pp_w_pp_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_u_pp_w_pp_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_u_pp_w_pp_avg_local[i] = 0.0;
            rho_u_pp_w_pp_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
                double* u = data_velocity->getPointer(0);
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double u_pp = u[idx_velocity] - rho_u_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    const double w_pp = w[idx_velocity] - rho_w_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    rho_u_pp_w_pp_avg_local[idx_fine] +=
                                        rho[idx_density]*u_pp*w_pp*weight_local;
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
        
        mpi.Reduce(
            rho_u_pp_w_pp_avg_local,
            rho_u_pp_w_pp_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean density multiplied by Reynolds shear stress component
         * in x- and z-directions inside mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double rho_R13_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    rho_R13_sum += rho_u_pp_w_pp_avg_global[i]*dx;
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_R13_sum;
            }
            else
            {
                double rho_R13_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        rho_R13_sum += rho_u_pp_w_pp_avg_global[i];
                        count++;
                    }
                }
                
                const double rho_R13_mean = rho_R13_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_R13_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_u_avg_local);
        std::free(rho_u_avg_global);
        std::free(rho_w_avg_local);
        std::free(rho_w_avg_global);
        std::free(rho_u_pp_w_pp_avg_local);
        std::free(rho_u_pp_w_pp_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean Reynolds shear stress component in y- and z-directions inside mixing layer
 * with assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputReynoldsShearStressInYZDirectionsInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'R23_INT_HOMO_YZ'/'R23_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'R23_INT_HOMO_YZ'/'R23_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'R23_INT_HOMO_YZ'/'R23_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
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
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_v_to_add = rho[idx_density]*v[idx_velocity]*weight_local;
                                const double rho_w_to_add = rho[idx_density]*w[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
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
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
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
        
        double* rho_v_pp_w_pp_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_pp_w_pp_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_v_pp_w_pp_avg_local[i] = 0.0;
            rho_v_pp_w_pp_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double v_pp = v[idx_velocity] - rho_v_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    const double w_pp = w[idx_velocity] - rho_w_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    rho_v_pp_w_pp_avg_local[idx_fine] +=
                                        rho[idx_density]*v_pp*w_pp*weight_local;
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
        
        mpi.Reduce(
            rho_v_pp_w_pp_avg_local,
            rho_v_pp_w_pp_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean Reynolds shear stress component in y- and z-directions
         * inside mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double R23_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    R23_sum += rho_v_pp_w_pp_avg_global[i]*dx/rho_avg_global[i];
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << R23_sum;
            }
            else
            {
                double R23_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        R23_sum += rho_v_pp_w_pp_avg_global[i]/rho_avg_global[i];
                        count++;
                    }
                }
                
                const double R23_mean = R23_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << R23_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_v_avg_local);
        std::free(rho_v_avg_global);
        std::free(rho_w_avg_local);
        std::free(rho_w_avg_global);
        std::free(rho_v_pp_w_pp_avg_local);
        std::free(rho_v_pp_w_pp_avg_global);
    }
    
    if (mpi.getRank() == 0)
    {
        f_out.close();
    }
}


/*
 * Output integrated or mean density multiplied by Reynolds shear stress component in y- and z-directions
 * inside mixing layer with assumed homogeneity in yz-plane to a file.
 */
void
SecondMomentStatisticsUtilities::
outputDensityTimesReynoldsShearStressInYZDirectionsInMixingLayerWithHomogeneityInYZPlane(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const bool integrated)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    if (d_num_species != 2)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "'rR23_INT_HOMO_YZ'/'rR23_HOMO_YZ_IN_ML_X' can be computed with two species only."
            << std::endl);
    }
    
    if (d_flow_model.expired())
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "The object is not setup yet!"
            << std::endl);
    }
    
    boost::shared_ptr<FlowModel> flow_model_tmp = d_flow_model.lock();
    
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
            << "There is no 'rR23_INT_HOMO_YZ'/'rR23_HOMO_YZ_IN_ML_X' for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "There is no 'rR23_INT_HOMO_YZ'/'rR23_HOMO_YZ_IN_ML_X' for two-dimensional problem."
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
        
        double* X_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* X_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_w_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            X_avg_local[i] = 0.0;
            X_avg_global[i] = 0.0;
            rho_avg_local[i] = 0.0;
            rho_avg_global[i] = 0.0;
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
                 * Register the patch, mole fractions, density and velocity in the flow model and compute
                 * the corresponding cell data.
                 */
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("MOLE_FRACTIONS", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to mole fractions, density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_mole_fractions =
                    flow_model_tmp->getCellData("MOLE_FRACTIONS");
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                std::vector<double*> X;
                X.reserve(d_num_species);
                for (int si = 0; si < d_num_species; si++)
                {
                    X.push_back(data_mole_fractions->getPointer(si));
                }
                double* rho = data_density->getPointer(0);
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
                
                const hier::IntVector num_ghosts_mole_fractions = data_mole_fractions->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_mole_fractions = data_mole_fractions->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_density = data_density->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_density = data_density->getGhostBox().numberCells();
                
                const hier::IntVector num_ghosts_velocity = data_velocity->getGhostCellWidth();
                const hier::IntVector ghostcell_dims_velocity = data_velocity->getGhostBox().numberCells();
                
                const int num_ghosts_0_mole_fractions = num_ghosts_mole_fractions[0];
                const int num_ghosts_1_mole_fractions = num_ghosts_mole_fractions[1];
                const int num_ghosts_2_mole_fractions = num_ghosts_mole_fractions[2];
                const int ghostcell_dim_0_mole_fractions = ghostcell_dims_mole_fractions[0];
                const int ghostcell_dim_1_mole_fractions = ghostcell_dims_mole_fractions[1];
                
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
                                
                                const int idx_mole_fractions = (relative_idx_lo_0 + i + num_ghosts_0_mole_fractions) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_mole_fractions)*ghostcell_dim_0_mole_fractions +
                                    (relative_idx_lo_2 + k + num_ghosts_2_mole_fractions)*ghostcell_dim_0_mole_fractions*
                                        ghostcell_dim_1_mole_fractions;
                                
                                const int idx_density = (relative_idx_lo_0 + i + num_ghosts_0_density) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_density)*ghostcell_dim_0_density +
                                    (relative_idx_lo_2 + k + num_ghosts_2_density)*ghostcell_dim_0_density*
                                        ghostcell_dim_1_density;
                                
                                const int idx_velocity = (relative_idx_lo_0 + i + num_ghosts_0_velocity) +
                                    (relative_idx_lo_1 + j + num_ghosts_1_velocity)*ghostcell_dim_0_velocity +
                                    (relative_idx_lo_2 + k + num_ghosts_2_velocity)*ghostcell_dim_0_velocity*
                                        ghostcell_dim_1_velocity;
                                
                                const double weight_local = weight/((double) n_overlapped);
                                
                                const double X_to_add = X[0][idx_mole_fractions]*weight_local;
                                const double rho_to_add = rho[idx_density]*weight_local;
                                const double rho_v_to_add = rho[idx_density]*v[idx_velocity]*weight_local;
                                const double rho_w_to_add = rho[idx_density]*w[idx_velocity]*weight_local;
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    X_avg_local[idx_fine] += X_to_add;
                                    rho_avg_local[idx_fine] += rho_to_add;
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
                
                flow_model_tmp->unregisterPatch();
            }
        }
        
        mpi.Allreduce(
            X_avg_local,
            X_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM);
        
        mpi.Allreduce(
            rho_avg_local,
            rho_avg_global,
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
        
        double* rho_v_pp_w_pp_avg_local = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        double* rho_v_pp_w_pp_avg_global = (double*)std::malloc(finest_level_dim_0*sizeof(double));
        
        for (int i = 0; i < finest_level_dim_0; i++)
        {
            rho_v_pp_w_pp_avg_local[i] = 0.0;
            rho_v_pp_w_pp_avg_global[i] = 0.0;
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
                
                flow_model_tmp->registerPatchWithDataContext(*patch, data_context);
                
                hier::IntVector num_ghosts = flow_model_tmp->getNumberOfGhostCells();
                
                std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("DENSITY", hier::IntVector::getZero(d_dim)));
                
                num_subghosts_of_data.insert(
                    std::pair<std::string, hier::IntVector>("VELOCITY", hier::IntVector::getZero(d_dim)));
                
                flow_model_tmp->registerDerivedVariables(num_subghosts_of_data);
                
                flow_model_tmp->computeDerivedCellData();
                
                /*
                 * Get the pointers to density and velocity data inside the flow model.
                 */
                
                boost::shared_ptr<pdat::CellData<double> > data_density =
                    flow_model_tmp->getCellData("DENSITY");
                
                boost::shared_ptr<pdat::CellData<double> > data_velocity =
                    flow_model_tmp->getCellData("VELOCITY");
                
                double* rho = data_density->getPointer(0);
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
                                
                                for (int ii = 0; ii < ratioToFinestLevel_0; ii++)
                                {
                                    const int idx_fine = (idx_lo_0 + i)*ratioToFinestLevel_0 + ii;
                                    
                                    const double v_pp = v[idx_velocity] - rho_v_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    const double w_pp = w[idx_velocity] - rho_w_avg_global[idx_fine]/
                                        rho_avg_global[idx_fine];
                                    
                                    rho_v_pp_w_pp_avg_local[idx_fine] +=
                                        rho[idx_density]*v_pp*w_pp*weight_local;
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
        
        mpi.Reduce(
            rho_v_pp_w_pp_avg_local,
            rho_v_pp_w_pp_avg_global,
            finest_level_dim_0,
            MPI_DOUBLE,
            MPI_SUM,
            0);
        
        /*
         * Compute and output the integrated or mean density multiplied by Reynolds shear stress component
         * in y- and z-directions inside mixing layer (only done by process 0).
         */
        
        if (mpi.getRank() == 0)
        {
            if (integrated == true)
            {
                const int finest_level_dim_0 = finest_level_dims[0];
                
                const double dx = (x_hi[0] - x_lo[0])/(double(finest_level_dim_0));
                
                double rho_R23_sum = 0.0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    rho_R23_sum += rho_v_pp_w_pp_avg_global[i]*dx;
                }
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_R23_sum;
            }
            else
            {
                double rho_R23_sum = 0.0;
                int count = 0;
                
                for (int i = 0; i < finest_level_dim_0; i++)
                {
                    if (4.0*X_avg_global[i]*(1.0 - X_avg_global[i]) > 0.9)
                    {
                        rho_R23_sum += rho_v_pp_w_pp_avg_global[i];
                        count++;
                    }
                }
                
                const double rho_R23_mean = rho_R23_sum/count;
                
                f_out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)
                      << "\t" << rho_R23_mean;
            }
        }
        
        std::free(X_avg_local);
        std::free(X_avg_global);
        std::free(rho_avg_local);
        std::free(rho_avg_global);
        std::free(rho_v_avg_local);
        std::free(rho_v_avg_global);
        std::free(rho_w_avg_local);
        std::free(rho_w_avg_global);
        std::free(rho_v_pp_w_pp_avg_local);
        std::free(rho_v_pp_w_pp_avg_global);
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
SecondMomentStatisticsUtilities::outputNumberOfCells(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
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
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
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
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
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
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
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
SecondMomentStatisticsUtilities::outputWeightedNumberOfCells(
    const std::string& stat_dump_filename,
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
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
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
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
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
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
            
            boost::shared_ptr<hier::PatchLevel> patch_level(
                patch_hierarchy->getPatchLevel(li));
            
            for (hier::PatchLevel::iterator ip(patch_level->begin());
                 ip != patch_level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch> patch = *ip;
                
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
            
            if (statistical_quantity_key == "MIXING_WIDTH_MOL_FR_X")
            {
                f_out << "\t" << "MIXING_WIDTH_MOL_FR_X";
            }
            else if (statistical_quantity_key == "MIXEDNESS_MOL_FR_X")
            {
                f_out << "\t" << "MIXEDNESS_MOL_FR_X   ";
            }
            else if (statistical_quantity_key == "a1_INT_HOMO_YZ")
            {
                f_out << "\t" << "a1_INT_HOMO_YZ       ";
            }
            else if (statistical_quantity_key == "a1_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "a1_HOMO_YZ_IN_ML_X   ";
            }
            else if (statistical_quantity_key == "ra1_INT_HOMO_YZ")
            {
                f_out << "\t" << "ra1_INT_HOMO_YZ      ";
            }
            else if (statistical_quantity_key == "ra1_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "ra1_HOMO_YZ_IN_ML_X  ";
            }
            else if (statistical_quantity_key == "b_INT_HOMO_YZ")
            {
                f_out << "\t" << "b_INT_HOMO_YZ        ";
            }
            else if (statistical_quantity_key == "b_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "b_HOMO_YZ_IN_ML_X    ";
            }
            else if (statistical_quantity_key == "rb_INT_HOMO_YZ")
            {
                f_out << "\t" << "rb_INT_HOMO_YZ       ";
            }
            else if (statistical_quantity_key == "rb_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "rb_HOMO_YZ_IN_ML_X   ";
            }
            else if (statistical_quantity_key == "R11_INT_HOMO_YZ")
            {
                f_out << "\t" << "R11_INT_HOMO_YZ      ";
            }
            else if (statistical_quantity_key == "R11_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "R11_HOMO_YZ_IN_ML_X  ";
            }
            else if (statistical_quantity_key == "rR11_INT_HOMO_YZ")
            {
                f_out << "\t" << "rR11_INT_HOMO_YZ     ";
            }
            else if (statistical_quantity_key == "rR11_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "rR11_HOMO_YZ_IN_ML_X ";
            }
            else if (statistical_quantity_key == "R22_INT_HOMO_YZ")
            {
                f_out << "\t" << "R22_INT_HOMO_YZ      ";
            }
            else if (statistical_quantity_key == "R22_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "R22_HOMO_YZ_IN_ML_X  ";
            }
            else if (statistical_quantity_key == "rR22_INT_HOMO_YZ")
            {
                f_out << "\t" << "rR22_INT_HOMO_YZ     ";
            }
            else if (statistical_quantity_key == "rR22_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "rR22_HOMO_YZ_IN_ML_X ";
            }
            else if (statistical_quantity_key == "R33_INT_HOMO_YZ")
            {
                f_out << "\t" << "R33_INT_HOMO_YZ      ";
            }
            else if (statistical_quantity_key == "R33_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "R33_HOMO_YZ_IN_ML_X  ";
            }
            else if (statistical_quantity_key == "rR33_INT_HOMO_YZ")
            {
                f_out << "\t" << "rR33_INT_HOMO_YZ     ";
            }
            else if (statistical_quantity_key == "rR33_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "rR33_HOMO_YZ_IN_ML_X ";
            }
            else if (statistical_quantity_key == "R12_INT_HOMO_YZ")
            {
                f_out << "\t" << "R12_INT_HOMO_YZ      ";
            }
            else if (statistical_quantity_key == "R12_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "R12_HOMO_YZ_IN_ML_X  ";
            }
            else if (statistical_quantity_key == "rR12_INT_HOMO_YZ")
            {
                f_out << "\t" << "rR12_INT_HOMO_YZ     ";
            }
            else if (statistical_quantity_key == "rR12_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "rR12_HOMO_YZ_IN_ML_X ";
            }
            else if (statistical_quantity_key == "R13_INT_HOMO_YZ")
            {
                f_out << "\t" << "R13_INT_HOMO_YZ      ";
            }
            else if (statistical_quantity_key == "R13_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "R13_HOMO_YZ_IN_ML_X  ";
            }
            else if (statistical_quantity_key == "rR13_INT_HOMO_YZ")
            {
                f_out << "\t" << "rR13_INT_HOMO_YZ     ";
            }
            else if (statistical_quantity_key == "rR13_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "rR13_HOMO_YZ_IN_ML_X ";
            }
            else if (statistical_quantity_key == "R23_INT_HOMO_YZ")
            {
                f_out << "\t" << "R23_INT_HOMO_YZ      ";
            }
            else if (statistical_quantity_key == "R23_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "R23_HOMO_YZ_IN_ML_X  ";
            }
            else if (statistical_quantity_key == "rR23_INT_HOMO_YZ")
            {
                f_out << "\t" << "rR23_INT_HOMO_YZ     ";
            }
            else if (statistical_quantity_key == "rR23_HOMO_YZ_IN_ML_X")
            {
                f_out << "\t" << "rR23_HOMO_YZ_IN_ML_X ";
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
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const boost::shared_ptr<hier::VariableContext>& data_context)
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!stat_dump_filename.empty());
#endif
    
    boost::shared_ptr<SecondMomentStatisticsUtilities> second_moment_statistics_utilities(
        new SecondMomentStatisticsUtilities(
            "second moment statistics utilities",
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
        
        if (statistical_quantity_key == "MIXING_WIDTH_MOL_FR_X")
        {
            second_moment_statistics_utilities->outputMixingWidthInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "MIXEDNESS_MOL_FR_X")
        {
            second_moment_statistics_utilities->outputMixednessInXDirection(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "a1_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputTurbulentMassFluxVelocityInXDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "a1_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputTurbulentMassFluxVelocityInXDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "ra1_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputTurbulentMassFluxInXDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "ra1_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputTurbulentMassFluxInXDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "b_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputDensitySpecificVolumeCovarianceInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "b_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputDensitySpecificVolumeCovarianceInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "rb_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputDensityTimesDensitySpecificVolumeCovarianceInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "rb_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputDensityTimesDensitySpecificVolumeCovarianceInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "R11_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputReynoldsNormalStressInXDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "R11_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputReynoldsNormalStressInXDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "rR11_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputDensityTimesReynoldsNormalStressInXDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "rR11_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputDensityTimesReynoldsNormalStressInXDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "R22_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputReynoldsNormalStressInYDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "R22_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputReynoldsNormalStressInYDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "rR22_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputDensityTimesReynoldsNormalStressInYDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "rR22_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputDensityTimesReynoldsNormalStressInYDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "R33_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputReynoldsNormalStressInZDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "R33_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputReynoldsNormalStressInZDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "rR33_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputDensityTimesReynoldsNormalStressInZDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "rR33_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputDensityTimesReynoldsNormalStressInZDirectionInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "R12_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputReynoldsShearStressInXYDirectionsInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "R12_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputReynoldsShearStressInXYDirectionsInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "rR12_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputDensityTimesReynoldsShearStressInXYDirectionsInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "rR12_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputDensityTimesReynoldsShearStressInXYDirectionsInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "R13_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputReynoldsShearStressInXZDirectionsInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "R13_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputReynoldsShearStressInXZDirectionsInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "rR13_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputDensityTimesReynoldsShearStressInXZDirectionsInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "rR13_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputDensityTimesReynoldsShearStressInXZDirectionsInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "R23_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputReynoldsShearStressInYZDirectionsInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "R23_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputReynoldsShearStressInYZDirectionsInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "rR23_INT_HOMO_YZ")
        {
            second_moment_statistics_utilities->
                outputDensityTimesReynoldsShearStressInYZDirectionsInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    true);
        }
        else if (statistical_quantity_key == "rR23_HOMO_YZ_IN_ML_X")
        {
            second_moment_statistics_utilities->
                outputDensityTimesReynoldsShearStressInYZDirectionsInMixingLayerWithHomogeneityInYZPlane(
                    stat_dump_filename,
                    patch_hierarchy,
                    data_context,
                    false);
        }
        else if (statistical_quantity_key == "NUM_CELLS")
        {
            second_moment_statistics_utilities->outputNumberOfCells(
                stat_dump_filename,
                patch_hierarchy,
                data_context);
        }
        else if (statistical_quantity_key == "WEIGHTED_NUM_CELLS")
        {
            second_moment_statistics_utilities->outputWeightedNumberOfCells(
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
}


