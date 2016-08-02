#ifndef FLOW_MODEL_SINGLE_SPECIES_HPP
#define FLOW_MODEL_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModel.hpp"

#include "flow/flow_models/Riemann_solvers/RiemannSolverSingleSpeciesHLLC.hpp"
#include "flow/flow_models/Riemann_solvers/RiemannSolverSingleSpeciesHLLC_HLL.hpp"

class FlowModelSingleSpecies: public FlowModel
{
    public:
        FlowModelSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const hier::IntVector& num_ghosts,
            const int& num_species,
            const boost::shared_ptr<EquationOfState>& equation_of_state);
        
        ~FlowModelSingleSpecies() {}
        
        /*
         * Print all characteristics of the flow model class.
         */
        void printClassData(std::ostream& os) const;
        
        /*
         * Register the conservative variables.
         */
        void
        registerConservativeVariables(RungeKuttaLevelIntegrator* integrator);
        
        /*
         * Get the names of conservative variables.
         */
        std::vector<std::string> getNamesOfConservativeVariables(bool have_underscores = false);
        
        /*
         * Get the names of primitive variables.
         */
        std::vector<std::string> getNamesOfPrimitiveVariables(bool have_underscores = false);
        
        /*
         * Get the variable types of conservative variables.
         */
        std::vector<std::string> getVariableTypesOfConservativeVariables();
        
        /*
         * Get the variable types of primitive variables.
         */
        std::vector<std::string> getVariableTypesOfPrimitiveVariables();
        
        /*
         * Get the conservative variables.
         */
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > >
        getConservativeVariables();
        
        /*
         * Register a patch with a data context.
         */
        void
        registerPatchWithDataContext(
            const hier::Patch& patch,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Register different derived variables in the registered patch. The derived variables to be registered
         * are given as entires in a map of the variable name to the number of sub-ghost cells required.
         * If the variable to be registered is one of the conservative variable, the corresponding entry
         * in the map is ignored.
         */
        void
        registerDerivedCellVariable(
            const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data);
        
        /*
         * Register the required variables for the computation of projection matrix
         * of conservative variables and its inverse at faces in the registered patch.
         */
        void
        registerFaceProjectionMatricesOfConservativeVariables(
            const hier::IntVector& num_subghosts,
            const AVERAGING& averaging);
        
        /*
         * Register the required variables for the computation of projection matrix
         * of primitive variables and its inverse at faces in the registered patch.
         */
        void
        registerFaceProjectionMatricesOfPrimitiveVariables(
            const hier::IntVector& num_subghosts,
            const AVERAGING& averaging);
        
        /*
         * Unregister the registered patch. The registered data context and all global derived
         * cell data in the patch are dumped.
         */
        void unregisterPatch();
        
        /*
         * Compute global cell data of different registered derived variables with the registered data context.
         */
        void
        computeGlobalDerivedCellData();
        
        /*
         * Get the global cell data of one cell variable in the registered patch.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getGlobalCellData(const std::string& variable_key);
        
        /*
         * Get the global cell data of different cell variables in the registered patch.
         */
        std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getGlobalCellData(const std::vector<std::string>& variable_keys);
        
        /*
         * Fill the interior global cell data of conservative variables with zeros.
         */
        void
        fillZeroGlobalCellDataConservativeVariables();
        
        /*
         * Update the interior global cell data of conservative variables.
         */
        void
        updateGlobalCellDataConservativeVariables();
        
        /*
         * Get the global cell data of the conservative variables in the registered patch.
         */
        std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getGlobalCellDataConservativeVariables();
        
        /*
         * Get the global cell data of the primitive variables in the registered patch.
         */
        std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getGlobalCellDataPrimitiveVariables();
        
        /*
         * Compute the local face data of projection matrix of conservative variables in the
         * registered patch.
         */
        void
        computeLocalFaceProjectionMatrixOfConservativeVariables(
            boost::multi_array<double, 2>& projection_matrix,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction);
        
        /*
         * Compute the local face data of inverse of projection matrix of conservative variables
         * in the registered patch.
         */
        void
        computeLocalFaceProjectionMatrixInverseOfConservativeVariables(
            boost::multi_array<double, 2>& projection_matrix_inv,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction);
        
        /*
         * Compute the local face data of projection matrix of primitive variables in the
         * registered patch.
         */
        void
        computeLocalFaceProjectionMatrixOfPrimitiveVariables(
            boost::multi_array<double, 2>& projection_matrix,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction);
        
        /*
         * Compute the local face data of inverse of projection matrix of primitive variables
         * in the registered patch.
         */
        void
        computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables(
            boost::multi_array<double, 2>& projection_matrix_inv,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction);
        
        /*
         * Compute the local intercell quantities with conservative variables on each side of the face
         * from Riemann solver at face.
         * fluxes_face: Convective flux at face.
         * velocity_face: Velocity at face.
         * The FlowModelSingleSpecies class modifies nothing for velocity_face.
         */
        void
        computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& flux_face,
            std::vector<boost::reference_wrapper<double> >& velocity_face,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_plus,
            const DIRECTION& direction,
            const RIEMANN_SOLVER& Riemann_solver);
        
        /*
         * Compute the local intercell quantities with primitive variables on each side of the face
         * from Riemann solver at face.
         * fluxes_face: Convective flux at face.
         * velocity_face: Velocity at face.
         * The FlowModelSingleSpecies class modifies nothing for velocity_face.
         */
        void
        computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& flux_face,
            std::vector<boost::reference_wrapper<double> >& velocity_face,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_plus,
            const DIRECTION& direction,
            const RIEMANN_SOLVER& Riemann_solver);
        
        /*
         * Check whether the given conservative variables are within the bounds.
         */
        bool
        haveConservativeVariablesBounded(const std::vector<double>& conservative_variables);
        
        /*
         * Check whether the given primitive variables are within the bounds.
         */
        bool
        havePrimitiveVariablesBounded(const std::vector<double>& primitive_variables);
        
        /*
         * Convert vector of pointers of conservative cell data to vectors of pointers of primitive cell data.
         */
        void
        convertLocalCellDataPointersConservativeVariablesToPrimitiveVariables(
            const std::vector<const double*>& conservative_variables,
            const std::vector<double*>& primitive_variables);
        
        /*
         * Convert vector of pointers of primitive cell data to vectors of pointers of conservative cell data.
         */
        void
        convertLocalCellDataPointersPrimitiveVariablesToConservativeVariables(
            const std::vector<const double*>& primitive_variables,
            const std::vector<double*>& conservative_variables);
        
        /*
         * Compute derived plot quantities registered with the VisIt data writers from data that
         * is maintained on each patch in the hierarchy.
         */
        bool
        packDerivedDataIntoDoubleBuffer(
            double* buffer,
            const hier::Patch& patch,
            const hier::Box& region,
            const std::string& variable_name,
            int depth_id,
            double simulation_time) const;
        
        /*
         * Register the plotting quantities.
         */
#ifdef HAVE_HDF5
        void
        registerPlotQuantities(
            const boost::shared_ptr<appu::VisItDataWriter>& visit_writer);
#endif
        
    private:
        /*
         * Set the ghost boxes and their dimensions of derived cell variables.
         */
        void
        setGhostBoxesAndDimensionsDerivedCellVariables();
        
        /*
         * Get the global cell data of density in the registered patch.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getGlobalCellDataDensity();
        
        /*
         * Get the global cell data of momentum in the registered patch.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getGlobalCellDataMomentum();
        
        /*
         * Get the global cell data of total energy in the registered patch.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getGlobalCellDataTotalEnergy();
        
        /*
         * Compute the global cell data of pressure in the registered patch.
         */
        void computeGlobalCellDataPressure();
        
        /*
         * Compute the global cell data of velocity in the registered patch.
         */
        void computeGlobalCellDataVelocity();
        
        /*
         * Compute the global cell data of sound speed with pressure in the registered patch.
         */
        void computeGlobalCellDataSoundSpeedWithPressure();
        
        /*
         * Compute the global cell data of dilatation with velocity in the registered patch.
         */
        void computeGlobalCellDataDilatationWithVelocity();
        
        /*
         * Compute the global cell data of vorticity with velocity in the registered patch.
         */
        void computeGlobalCellDataVorticityWithVelocity();
        
        /*
         * Compute the global cell data of enstrophy with velocity and vorticity in the registered patch.
         */
        void computeGlobalCellDataEnstrophyWithVelocityAndVorticity();
        
        /*
         * Compute the global cell data of convective flux with pressure and velocity in the registered patch.
         */
        void computeGlobalCellDataConvectiveFluxWithPressureAndVelocity(DIRECTION direction);
        
        /*
         * Compute the global cell data of maximum wave speed with pressure, velocity and sound speed in the registered patch.
         */
        void computeGlobalCellDataMaxWaveSpeedWithPressureVelocityAndSoundSpeed(DIRECTION direction);
        
        /*
         * boost::shared_ptr to registered conservative variables.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_momentum;
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_total_energy;
        
        /*
         * Number of sub-ghost cells of derived cell data.
         */
        hier::IntVector d_num_subghosts_pressure;
        hier::IntVector d_num_subghosts_velocity;
        hier::IntVector d_num_subghosts_sound_speed;
        hier::IntVector d_num_subghosts_dilatation;
        hier::IntVector d_num_subghosts_vorticity;
        hier::IntVector d_num_subghosts_enstrophy;
        hier::IntVector d_num_subghosts_convective_flux_x;
        hier::IntVector d_num_subghosts_convective_flux_y;
        hier::IntVector d_num_subghosts_convective_flux_z;
        hier::IntVector d_num_subghosts_max_wave_speed_x;
        hier::IntVector d_num_subghosts_max_wave_speed_y;
        hier::IntVector d_num_subghosts_max_wave_speed_z;
        
        /*
         * Boxes with sub-ghost cells of derived cell data.
         */
        hier::Box d_subghost_box_pressure;
        hier::Box d_subghost_box_velocity;
        hier::Box d_subghost_box_sound_speed;
        hier::Box d_subghost_box_dilatation;
        hier::Box d_subghost_box_vorticity;
        hier::Box d_subghost_box_enstrophy;
        hier::Box d_subghost_box_convective_flux_x;
        hier::Box d_subghost_box_convective_flux_y;
        hier::Box d_subghost_box_convective_flux_z;
        hier::Box d_subghost_box_max_wave_speed_x;
        hier::Box d_subghost_box_max_wave_speed_y;
        hier::Box d_subghost_box_max_wave_speed_z;
        
        /*
         * Dimensions of boxes with sub-ghost cells of derived cell data.
         */
        hier::IntVector d_subghostcell_dims_pressure;
        hier::IntVector d_subghostcell_dims_velocity;
        hier::IntVector d_subghostcell_dims_sound_speed;
        hier::IntVector d_subghostcell_dims_dilatation;
        hier::IntVector d_subghostcell_dims_vorticity;
        hier::IntVector d_subghostcell_dims_enstrophy;
        hier::IntVector d_subghostcell_dims_convective_flux_x;
        hier::IntVector d_subghostcell_dims_convective_flux_y;
        hier::IntVector d_subghostcell_dims_convective_flux_z;
        hier::IntVector d_subghostcell_dims_max_wave_speed_x;
        hier::IntVector d_subghostcell_dims_max_wave_speed_y;
        hier::IntVector d_subghostcell_dims_max_wave_speed_z;
        
        /*
         * boost::shared_ptr to derived cell data.
         */
        boost::shared_ptr<pdat::CellData<double> > d_data_pressure;
        boost::shared_ptr<pdat::CellData<double> > d_data_velocity;
        boost::shared_ptr<pdat::CellData<double> > d_data_sound_speed;
        boost::shared_ptr<pdat::CellData<double> > d_data_dilatation;
        boost::shared_ptr<pdat::CellData<double> > d_data_vorticity;
        boost::shared_ptr<pdat::CellData<double> > d_data_enstrophy;
        boost::shared_ptr<pdat::CellData<double> > d_data_convective_flux_x;
        boost::shared_ptr<pdat::CellData<double> > d_data_convective_flux_y;
        boost::shared_ptr<pdat::CellData<double> > d_data_convective_flux_z;
        boost::shared_ptr<pdat::CellData<double> > d_data_max_wave_speed_x;
        boost::shared_ptr<pdat::CellData<double> > d_data_max_wave_speed_y;
        boost::shared_ptr<pdat::CellData<double> > d_data_max_wave_speed_z;
        
        /*
         * Riemann solvers.
         */
        RiemannSolverSingleSpeciesHLLC     d_Riemann_solver_HLLC;
        RiemannSolverSingleSpeciesHLLC_HLL d_Riemann_solver_HLLC_HLL;
        
};

#endif /* FLOW_MODEL_SINGLE_SPECIES_HPP */
