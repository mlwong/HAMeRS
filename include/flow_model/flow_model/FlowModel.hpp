#ifndef FLOW_MODEL_HPP
#define FLOW_MODEL_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceData.h"

#include "equation_of_state/EquationOfStateIdealGas.hpp"
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

class FlowModel
{
    public:
        FlowModel(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const hier::IntVector& num_ghosts,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<EquationOfState>& equation_of_state):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_ghosts(num_ghosts),
                d_num_eqn(num_eqn),
                d_num_species(num_species),
                d_equation_of_state(equation_of_state),
                d_patch(nullptr),
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
         * Register a patch with global cell data of different variables in the patch.
         */
        virtual void
        registerPatchWithGlobalCellData(
            const hier::Patch& patch,
            const std::unordered_map<std::string, hier::IntVector>& num_subghosts_of_data,
            const boost::shared_ptr<hier::VariableContext>& data_context) = 0;
        
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
         * Unregister the registered patch and all global cell data in the patch.
         */
        virtual void unregisterPatchWithGlobalCellData() = 0;
        
        /*
         * Compute the global cell data of the registered variables in the registered patch.
         */
        virtual void
        computeGlobalCellData() = 0;
        
        /*
         * Get the global cell data of one cell variable in the registered patch.
         * The number of sub-ghost cells and the dimensions of box with sub-ghost cells are also returned.
         */
        virtual boost::shared_ptr<pdat::CellData<double> >
        getGlobalCellData(
            const std::string& variable_key,
            hier::IntVector& num_subghosts,
            hier::IntVector& subghostcell_dims) = 0;
        
        /*
         * Get the global cell data of different cell variables in the registered patch.
         * The numbers of sub-ghost cells and the dimensions of boxes with sub-ghost cells are also returned.
         */
        virtual std::vector<boost::shared_ptr<pdat::CellData<double> > >
        getGlobalCellData(
            const std::vector<std::string>& variable_keys,
            std::vector<hier::IntVector>& num_subghosts,
            std::vector<hier::IntVector>& subghostcell_dims) = 0;
        
        /*
         * Get the pointers to the global cell data of the conservative variables in the registered patch.
         * The numbers of sub-ghost cells and the dimensions of boxes with sub-ghost cells are also returned.
         */
        virtual std::vector<double*>
        getGlobalCellDataPointerConservativeVariables(
            std::vector<hier::IntVector>& num_subghosts,
            std::vector<hier::IntVector>& subghostcell_dims) = 0;
        
        /*
         * Get the pointers to the global cell data of the primitive variables in the registered patch.
         * The numbers of sub-ghost cells and the dimensions of boxes with sub-ghost cells are also returned.
         */
        virtual std::vector<double*>
        getGlobalCellDataPointerPrimitiveVariables(
            std::vector<hier::IntVector>& num_subghosts,
            std::vector<hier::IntVector>& subghostcell_dims) = 0;
        
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
         * Set the cell variables if single-species flow model is chosen.
         */
        void
        setVariablesForSingleSpecies(
            const boost::shared_ptr<pdat::CellVariable<double> >& density,
            const boost::shared_ptr<pdat::CellVariable<double> >& momentum,
            const boost::shared_ptr<pdat::CellVariable<double> >& total_energy)
        {
            d_variable_density = density;
            d_variable_momentum = momentum;
            d_variable_total_energy = total_energy;
        }
        
        /*
         * Set the cell variables if four-equation multi-species conservative flow model
         * is chosen.
         */
        void
        setVariablesForFourEqnConservative(
            const boost::shared_ptr<pdat::CellVariable<double> >& partial_density,
            const boost::shared_ptr<pdat::CellVariable<double> >& momentum,
            const boost::shared_ptr<pdat::CellVariable<double> >& total_energy)
        {
            d_variable_partial_density = partial_density;
            d_variable_momentum = momentum;
            d_variable_total_energy = total_energy;
        }
        
        /*
         * Set the cell variables if four-equation multi-species flow model
         * by Shyue is chosen.
         */
        void
        setVariablesForFourEqnShyue(
            const boost::shared_ptr<pdat::CellVariable<double> >& density,
            const boost::shared_ptr<pdat::CellVariable<double> >& momentum,
            const boost::shared_ptr<pdat::CellVariable<double> >& total_energy,
            const boost::shared_ptr<pdat::CellVariable<double> >& mass_fraction)
        {
            d_variable_density = density;
            d_variable_momentum = momentum;
            d_variable_total_energy = total_energy;
            d_variable_mass_fraction = mass_fraction;
        }
        
        /*
         * Set the cell variables if five-equation multi-species flow model
         * by Allaire is chosen.
         */
        void
        setVariablesForFiveEqnAllaire(
            const boost::shared_ptr<pdat::CellVariable<double> >& partial_density,
            const boost::shared_ptr<pdat::CellVariable<double> >& momentum,
            const boost::shared_ptr<pdat::CellVariable<double> >& total_energy,
            const boost::shared_ptr<pdat::CellVariable<double> >& volume_fraction)
        {
            d_variable_partial_density = partial_density;
            d_variable_momentum = momentum;
            d_variable_total_energy = total_energy;
            d_variable_volume_fraction = volume_fraction;
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
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * boost::shared_ptr to EquationOfState.
         */
        const boost::shared_ptr<EquationOfState> d_equation_of_state;
        
        /*
         * Form of each equation.
         */
        std::vector<EQUATION_FORM> d_eqn_form;
        
        /*
         * pointer to registered patch.
         */
        const hier::Patch* d_patch;
        
        /*
         * boost::shared_ptr to patch data context.
         */
        boost::shared_ptr<hier::VariableContext> d_data_context;
        
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
         * boost::shared_ptr to registered cell variables.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_partial_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_momentum;
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_total_energy;
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_mass_fraction;
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_volume_fraction;
        
};

#endif /* FLOW_MODEL_HPP */
