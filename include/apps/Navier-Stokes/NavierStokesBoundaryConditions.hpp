#ifndef NAVIER_STOKES_BOUNDARY_CONDITIONS_HPP
#define NAVIER_STOKES_BOUNDARY_CONDITIONS_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModels.hpp"
#include "apps/Navier-Stokes/NavierStokesSpecialBoundaryConditions.hpp"
#include "util/basic_boundary_conditions/BasicBoundaryConditions.hpp"
#include "util/basic_boundary_conditions/BoundaryUtilityStrategy.hpp"
#include "util/basic_boundary_conditions/CartesianBoundaryDefines.hpp"

#include <map>
#include <string>
#include <vector>

using namespace SAMRAI;

class NavierStokesBoundaryConditions:
    public BoundaryUtilityStrategy
{
    public:
        NavierStokesBoundaryConditions(
            const std::string& object_name,
            const std::string& project_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const FLOW_MODEL::TYPE& flow_model_type,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& boundary_conditions_db,
            const bool& is_from_restart);
        
        ~NavierStokesBoundaryConditions() {};
        
        /*
         * Print all characteristics of the boundary conditions class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the boundary conditions class into the restart
         * database.
         */
        void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * This routine is a concrete implementation of the virtual function
         * in the base class BoundaryUtilityStrategy.  It reads DIRICHLET
         * boundary state values from the given database with the
         * given name string idenifier.  The integer location index
         * indicates the face (in 3D) or edge (in 2D) to which the boundary
         * condition applies.
         */
        void
        readDirichletBoundaryDataEntry(
            const HAMERS_SHARED_PTR<tbox::Database>& db,
            std::string& db_name,
            int bdry_location_index);
        
        /*
         * This routine is a concrete implementation of the virtual function
         * in the base class BoundaryUtilityStrategy.  It is a blank implementation
         * for the purposes of this class.
         */
        void
        readNeumannBoundaryDataEntry(
            const HAMERS_SHARED_PTR<tbox::Database>& db,
            std::string& db_name,
            int bdry_location_index);
        
        /*
         * Set the data in ghost cells corresponding to physical boundary
         * conditions.  Specific boundary conditions are determined by
         * information specified in input file and numerical routines.
         */
        void
        setPhysicalBoundaryConditions(
            hier::Patch& patch,
            const double fill_time,
            const hier::IntVector& ghost_width_to_fill,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
    private:
        std::vector<double>
        readPrimitiveDataEntry(
            HAMERS_SHARED_PTR<tbox::Database> db,
            const std::string& db_name);
        
        /*
         * Set defaults for boundary conditions. Set to bogus values
         * for error checking.
         */
        void
        setDefaultBoundaryConditions();
        
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Name of the project.
         */
        std::string d_project_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Flow model type.
         */
        const FLOW_MODEL::TYPE d_flow_model_type;
        
        /*
         * Flow model.
         */
        const HAMERS_SHARED_PTR<FlowModel> d_flow_model;
        
        /*
         * Boundary condition cases and boundary values.
         * Options are: FLOW, REFLECT, DIRICHLET and variants for nodes and edges.
         *
         * Input file values are read into these arrays.
         */
        std::vector<int> d_master_bdry_node_conds;
        std::vector<int> d_master_bdry_edge_conds; // Used in 2D and 3D only.
        std::vector<int> d_master_bdry_face_conds; // Used in 3D only.
        
        /*
         * Boundary condition cases for scalar and vector (i.e., depth > 1) variables.
         * These are post-processed input values and are passed to the boundary routines.
         */
        std::vector<int> d_scalar_bdry_node_conds;
        std::vector<int> d_vector_bdry_node_conds;
        
        std::vector<int> d_scalar_bdry_edge_conds; // Used in 2D and 3D only.
        std::vector<int> d_vector_bdry_edge_conds; // Used in 2D and 3D only.
        
        std::vector<int> d_scalar_bdry_face_conds; // Used in 3D only.
        std::vector<int> d_vector_bdry_face_conds; // Used in 3D only.
        
        std::vector<int> d_node_bdry_edge; // Used in 2D only.
        std::vector<int> d_edge_bdry_face; // Used in 3D only.
        std::vector<int> d_node_bdry_face; // Used in 3D only.
        
        /*
         * Vectors of node (1D), edge (2D) or face (3D) boundary values for DIRICHLET case.
         */
        std::vector<std::vector<double> > d_bdry_node_conservative_var;
        std::vector<std::vector<double> > d_bdry_edge_conservative_var;
        std::vector<std::vector<double> > d_bdry_face_conservative_var;
        
        /*
         * HAMERS_SHARED_PTR to the special boundary conditions.
         */
        HAMERS_SHARED_PTR<NavierStokesSpecialBoundaryConditions> d_Navier_Stokes_special_boundary_conditions;
        
};

#endif /* NAVIER_STOKES_BOUNDARY_CONDITIONS_HPP */
