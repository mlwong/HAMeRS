#ifndef FLOW_MODEL_HPP
#define FLOW_MODEL_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceData.h"

#include "algs/integrator/RungeKuttaLevelIntegrator.hpp"
#include "util/equations_of_state/EquationsOfState.hpp"
#include "Directions.hpp"

#include "boost/multi_array.hpp"
#include "boost/shared_ptr.hpp"
#include <string>
#include <unordered_map>
#include <vector>

using namespace SAMRAI;

enum EQUATION_FORM { CONSERVATIVE_EQN,
                     ADVECTIVE_EQN };

enum VARIABLE_TYPE { CONSERVATIVE_VAR,
                     PRIMITIVE_VAR };

enum AVERAGING { SIMPLE_AVG,
                 ROE_AVG };

enum RIEMANN_SOLVER { HLLC_RIEMANN_SOLVER,
                      HLLC_HLL_RIEMANN_SOLVER };

/*
 * The class should at least be able to register, compute and get the following cell data of the single-species/mixture:
 * DENSITY, MOMENTUM, TOTAL_ENERGY, PRESSURE, VELOCITY, SOUND_SPEED, DILATATION, VORTICITY, ENSTROPHY, PRIMITIVE_VARIABLES
 * CONVECTIVE_FLUX_X, CONVECTIVE_FLUX_Y, CONVECTIVE_FLUX_Z, MAX_WAVE_SPEED_X, MAX_WAVE_SPEED_Y, MAX_WAVE_SPEED_Z
 */
class FlowModel:
    public appu::VisDerivedDataStrategy
{
    public:
        FlowModel(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const hier::IntVector& num_ghosts,
            const int& num_species,
            const int& num_eqn,
            const boost::shared_ptr<EquationOfState>& equation_of_state):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_ghosts(num_ghosts),
                d_num_species(num_species),
                d_num_eqn(num_eqn),
                d_equation_of_state(equation_of_state),
                d_patch(nullptr),
                d_patch_registered(false),
                d_interior_box(hier::Box::getEmptyBox(dim)),
                d_ghost_box(hier::Box::getEmptyBox(dim)),
                d_interior_dims(hier::IntVector::getZero(d_dim)),
                d_ghostcell_dims(hier::IntVector::getZero(d_dim)),
                d_proj_mat_conservative_var_registered(false),
                d_proj_mat_primitive_var_registered(false),
                d_proj_mat_conservative_var_averaging(SIMPLE_AVG),
                d_proj_mat_primitive_var_averaging(SIMPLE_AVG)
        {}
        
        /*
         * Get the total number of equations.
         */
        const int& getNumberOfEquations() const
        {
            return d_num_eqn;
        }
        
        /*
         * Return the form of each equation.
         */
        const std::vector<EQUATION_FORM>& getEquationsForm() const
        {
            return d_eqn_form;
        }
        
        /*
         * Print all characteristics of the flow model class.
         */
        virtual void printClassData(std::ostream& os) const = 0;
        
        /*
         * Register the conservative variables.
         */
        virtual void
        registerConservativeVariables(RungeKuttaLevelIntegrator* integrator) = 0;
        
        /*
         * Get the names of conservative variables.
         */
        virtual std::vector<std::string> getNamesOfConservativeVariables(bool have_underscores = false) = 0;
        
        /*
         * Get the names of primitive variables.
         */
        virtual std::vector<std::string> getNamesOfPrimitiveVariables(bool have_underscores = false) = 0;
        
        /*
         * Get the variable types of conservative variables.
         */
        virtual std::vector<std::string> getVariableTypesOfConservativeVariables() = 0;
        
        /*
         * Get the variable types of primitive variables.
         */
        virtual std::vector<std::string> getVariableTypesOfPrimitiveVariables() = 0;
        
        /*
         * Get the conservative variables.
         */
        virtual std::vector<boost::shared_ptr<pdat::CellVariable<double> > >
        getConservativeVariables() = 0;
        
        /*
         * Register a patch with a data context.
         */
        virtual void
        registerPatchWithDataContext(
            const hier::Patch& patch,
            const boost::shared_ptr<hier::VariableContext>& data_context) = 0;
        
        /*
         * Register different derived variables in the registered patch. The derived variables to be registered
         * are given as entires in a map of the variable name to the number of sub-ghost cells required.
         * If the variable to be registered is one of the conservative variable, the corresponding entry
         * in the map is ignored.
         */
        virtual void
        registerDerivedCellVariable(
            const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data) = 0;
        
        /*
         * Register the required derived variables for the computation of projection matrix
         * of conservative variables and its inverse at faces in the registered patch.
         */
        virtual void
        registerFaceProjectionMatricesOfConservativeVariables(
            const hier::IntVector& num_subghosts,
            const AVERAGING& averaging) = 0;
        
        /*
         * Register the required derived variables for the computation of projection matrix
         * of primitive variables and its inverse at faces in the registered patch.
         */
        virtual void
        registerFaceProjectionMatricesOfPrimitiveVariables(
            const hier::IntVector& num_subghosts,
            const AVERAGING& averaging) = 0;
        
        /*
         * Unregister the registered patch. The registered data context and all global derived
         * cell data in the patch are dumped.
         */
        virtual void unregisterPatch() = 0;
        
        /*
         * Compute global cell data of different registered derived variables with the registered data context.
         */
        virtual void
        computeGlobalDerivedCellData() = 0;
        
        /*
         * Get the global cell data of one cell variable in the registered patch.
         */
        virtual boost::shared_ptr<pdat::CellData<double> >
        getGlobalCellData(const std::string& variable_key) = 0;
        
        /*
         * Get the global cell data of different cell variables in the registered patch.
         */
        virtual std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getGlobalCellData(const std::vector<std::string>& variable_keys) = 0;
        
        /*
         * Fill the interior global cell data of conservative variables with zeros.
         */
        virtual void
        fillZeroGlobalCellDataConservativeVariables() = 0;
        
        /*
         * Update the interior global cell data of conservative variables after time advancement.
         */
        virtual void
        updateGlobalCellDataConservativeVariables() = 0;
        
        /*
         * Get the global cell data of the conservative variables in the registered patch.
         */
        virtual std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getGlobalCellDataConservativeVariables() = 0;
        
        /*
         * Get the global cell data of the primitive variables in the registered patch.
         */
        virtual std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getGlobalCellDataPrimitiveVariables() = 0;
        
        /*
         * Compute the local face datum of projection matrix of conservative variables in the
         * registered patch.
         */
        virtual void
        computeLocalFaceProjectionMatrixOfConservativeVariables(
            boost::multi_array<double, 2>& projection_matrix,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction) = 0;
        
        /*
         * Compute the local face datum of inverse of projection matrix of conservative variables
         * in the registered patch.
         */
        virtual void
        computeLocalFaceProjectionMatrixInverseOfConservativeVariables(
            boost::multi_array<double, 2>& projection_matrix_inv,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction) = 0;
        
        /*
         * Compute the local face datum of projection matrix of primitive variables in the
         * registered patch.
         */
        virtual void
        computeLocalFaceProjectionMatrixOfPrimitiveVariables(
            boost::multi_array<double, 2>& projection_matrix,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction) = 0;
        
        /*
         * Compute the local face datum of inverse of projection matrix of primitive variables
         * in the registered patch.
         */
        virtual void
        computeLocalFaceProjectionMatrixInverseOfPrimitiveVariables(
            boost::multi_array<double, 2>& projection_matrix_inv,
            const hier::Index& cell_index_minus,
            const hier::Index& cell_index_plus,
            const DIRECTION& direction) = 0;
        
        /*
         * Compute the local intercell quantities with conservative variables on each side of the face
         * from Riemann solver at face.
         * fluxes_face: Convective flux at face.
         * velocity_face: Velocity at face.
         * The FlowModelSingleSpecies class modifies nothing for velocity_face.
         */
        virtual void
        computeLocalFaceFluxAndVelocityFromRiemannSolverWithConservativeVariables(
            std::vector<boost::reference_wrapper<double> >& flux_face,
            std::vector<boost::reference_wrapper<double> >& velocity_face,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& conservative_variables_plus,
            const DIRECTION& direction,
            const RIEMANN_SOLVER& Riemann_solver) = 0;
        
        /*
         * Compute the local intercell quantities with primitive variables on each side of the face
         * from Riemann solver at face.
         * fluxes_face: Convective flux at face.
         * velocity_face: Velocity at face.
         * The FlowModelSingleSpecies class modify nothing for velocity_face.
         */
        virtual void
        computeLocalFaceFluxAndVelocityFromRiemannSolverWithPrimitiveVariables(
            std::vector<boost::reference_wrapper<double> >& flux_face,
            std::vector<boost::reference_wrapper<double> >& velocity_face,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_minus,
            const std::vector<boost::reference_wrapper<double> >& primitive_variables_plus,
            const DIRECTION& direction,
            const RIEMANN_SOLVER& Riemann_solver) = 0;
        
        /*
         * Check whether the given conservative variables are within the bounds.
         */
        virtual bool
        haveConservativeVariablesBounded(const std::vector<double>& conservative_variables) = 0;
        
        /*
         * Check whether the given primitive variables are within the bounds.
         */
        virtual bool
        havePrimitiveVariablesBounded(const std::vector<double>& primitive_variables) = 0;
        
        /*
         * Convert vector of pointers of conservative cell data to vectors of pointers of primitive cell data.
         */
        virtual void
        convertLocalCellDataPointersConservativeVariablesToPrimitiveVariables(
            const std::vector<const double*>& conservative_variables,
            const std::vector<double*>& primitive_variables) = 0;
        
        /*
         * Convert vector of pointers of primitive cell data to vectors of pointers of conservative cell data.
         */
        virtual void
        convertLocalCellDataPointersPrimitiveVariablesToConservativeVariables(
            const std::vector<const double*>& primitive_variables,
            const std::vector<double*>& conservative_variables) = 0;
        
        /*
         * Set the plotting context.
         */
        void
        setPlotContext(
            const boost::shared_ptr<hier::VariableContext>& plot_context)
        {
            d_plot_context = plot_context;
        }
        
        /*
         * Compute derived plot quantities registered with the VisIt data writers from data that
         * is maintained on each patch in the hierarchy.
         */
        virtual bool
        packDerivedDataIntoDoubleBuffer(
            double* buffer,
            const hier::Patch& patch,
            const hier::Box& region,
            const std::string& variable_name,
            int depth_id,
            double simulation_time) const = 0;
        
        /*
         * Register the plotting quantities.
         */
#ifdef HAVE_HDF5
        virtual void
        registerPlotQuantities(
            const boost::shared_ptr<appu::VisItDataWriter>& visit_writer) = 0;
#endif
        
        /*
         * Set the equation of state used.
         */
        void
        setEquationOfState(const boost::shared_ptr<EquationOfState>& equation_of_state)
        {
            d_equation_of_state = equation_of_state;
        }
        
        /*
         * Set the number of ghost cells needed.
         */
        void
        setNumberOfGhostCells(const hier::IntVector& num_ghosts)
        {
            d_num_ghosts = num_ghosts;
        }
        
    protected:
        /*
         * Set the context for data on a patch.
         */
        void
        setDataContext(const boost::shared_ptr<hier::VariableContext>& context)
        {
           d_data_context = context;
        }
        
        /*
         * Return boost::shared_ptr to patch data context.
         */
        boost::shared_ptr<hier::VariableContext>
        getDataContext() const
        {
           return d_data_context;
        }
        
        /*
         * Reset the data context to be null.
         */
        void
        clearDataContext()
        {
           d_data_context.reset();
        }
        
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
         * Number of ghost cells for registered variables.
         */
        hier::IntVector d_num_ghosts;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * boost::shared_ptr to EquationOfState.
         */
        boost::shared_ptr<EquationOfState> d_equation_of_state;
        
        /*
         * Form of each equation.
         */
        std::vector<EQUATION_FORM> d_eqn_form;
        
        /*
         * Pointer to registered patch.
         */
        const hier::Patch* d_patch;
        
        /*
         * boost::shared_ptr to patch data context.
         */
        boost::shared_ptr<hier::VariableContext> d_data_context;
        
        /*
         * Boolean to determine whether a patch is already registered.
         */
        bool d_patch_registered;
        
        /*
         * Interior box and box with ghost cells.
         */
        hier::Box d_interior_box;
        hier::Box d_ghost_box;
        
        /*
         * Dimensions of interior box and box with ghost cells.
         */
        hier::IntVector d_interior_dims;
        hier::IntVector d_ghostcell_dims;
        
        /*
         * Properties of the projection matrix and its inverse.
         */
        bool d_proj_mat_conservative_var_registered;
        bool d_proj_mat_primitive_var_registered;
        AVERAGING d_proj_mat_conservative_var_averaging;
        AVERAGING d_proj_mat_primitive_var_averaging;
        
        /*
         * boost::shared_ptr to the plotting context.
         */
        boost::shared_ptr<hier::VariableContext> d_plot_context;
        
};

#endif /* FLOW_MODEL_HPP */
