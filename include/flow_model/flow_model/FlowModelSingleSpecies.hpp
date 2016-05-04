#ifndef FLOW_MODEL_SINGLE_SPECIES_HPP
#define FLOW_MODEL_SINGLE_SPECIES_HPP

#include "flow_model/flow_model/FlowModel.hpp"

#include "flow_model/Riemann_solver/RiemannSolverSingleSpeciesHLLC.hpp"
#include "flow_model/Riemann_solver/RiemannSolverSingleSpeciesHLLC_HLL.hpp"

class FlowModelSingleSpecies: public FlowModel
{
    public:
        FlowModelSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const hier::IntVector& num_ghosts,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<EquationOfState>& equation_of_state);
        
        /*
         * Print all characteristics of the flow model class.
         */
        void printClassData(std::ostream& os) const;
        
        /*
         * Register a patch with global cell data of different variables in the patch.
         */
        void
        registerPatchWithGlobalCellData(
            const hier::Patch& patch,
            const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Register the required derived variables for the computation of projection matrix
         * of conservative variables and its inverse at faces in the registered patch.
         */
        void
        registerFaceProjectionMatricesOfConservativeVariables(
            const hier::IntVector& num_subghosts,
            const AVERAGING& averaging);
        
        /*
         * Register the required derived variables for the computation of projection matrix
         * of primitive variables and its inverse at faces in the registered patch.
         */
        void
        registerFaceProjectionMatricesOfPrimitiveVariables(
            const hier::IntVector& num_subghosts,
            const AVERAGING& averaging);
        
        /*
         * Unregister the registered patch and all global cell data in the patch.
         */
        void unregisterPatchWithGlobalCellData();
        
        /*
         * Compute the global cell data of the registered variables in the registered patch.
         */
        void
        computeGlobalCellData();
        
        /*
         * Get the global cell data of one cell variable in the registered patch.
         * The number of sub-ghost cells and the dimensions of box with sub-ghost cells are also returned.
         */
        boost::shared_ptr<pdat::CellData<double> >
        getGlobalCellData(
            const std::string& variable_key,
            hier::IntVector& num_subghosts,
            hier::IntVector& subghostcell_dims);
        
        /*
         * Get the global cell data of different cell variables in the registered patch.
         * The numbers of sub-ghost cells and the dimensions of boxes with sub-ghost cells are also returned.
         */
        std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getGlobalCellData(
            const std::vector<std::string>& variable_keys,
            std::vector<hier::IntVector>& num_subghosts,
            std::vector<hier::IntVector>& subghostcell_dims);
        
        /*
         * Get the pointers to the global cell data of the conservative variables in the registered patch.
         * The numbers of sub-ghost cells and the dimensions of boxes with sub-ghost cells are also returned.
         */
        std::vector<double*>
        getGlobalCellDataPointerConservativeVariables(
            std::vector<hier::IntVector>& num_subghosts,
            std::vector<hier::IntVector>& subghostcell_dims);
        
        /*
         * Get the pointers to the global cell data of the primitive variables in the registered patch.
         * The numbers of sub-ghost cells and the dimensions of boxes with sub-ghost cells are also returned.
         */
        std::vector<double*>
        getGlobalCellDataPointerPrimitiveVariables(
            std::vector<hier::IntVector>& num_subghosts,
            std::vector<hier::IntVector>& subghostcell_dims);
        
        /*
         * Compute the local face datum of projection matrix of conservative variables in the
         * registered patch.
         */
        void
        computeLocalFaceProjectionMatrixOfConservativeVariables(
            boost::multi_array<double, 2>& projection_matrix,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction);
        
        /*
         * Compute the local face datum of inverse of projection matrix of conservative variables
         * in the registered patch.
         */
        void
        computeLocalFaceProjectionMatrixInverseOfConservativeVariables(
            boost::multi_array<double, 2>& projection_matrix_inv,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction);
        
        /*
         * Compute the local face datum of projection matrix of primitive variables in the
         * registered patch.
         */
        void
        computeLocalFaceProjectionMatrixOfPrimitiveVariables(
            boost::multi_array<double, 2>& projection_matrix,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction);
        
        /*
         * Compute the local face datum of inverse of projection matrix of primitive variables
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
        
    private:
        /*
         * Set the ghost boxes and their dimensions.
         */
        void
        setGhostBoxesAndDimensions();
        
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
        
        /*
         * Riemann solvers.
         */
        RiemannSolverSingleSpeciesHLLC     d_Riemann_solver_HLLC;
        RiemannSolverSingleSpeciesHLLC_HLL d_Riemann_solver_HLLC_HLL;
        
};

#endif /* FLOW_MODEL_SINGLE_SPECIES_HPP */
