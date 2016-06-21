#ifndef EULER_BOUNDARY_CONDITIONS_HPP
#define EULER_BOUNDARY_CONDITIONS_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

#include "flow/flow_models/FlowModels.hpp"

#include "boost/shared_ptr.hpp"
#include <map>
#include <string>
#include <vector>

using namespace SAMRAI;

/*
 * Function to print out enum NodeBdyLoc2D::Type value as text.
 */
inline std::ostream& operator<<(std::ostream& os, const NodeBdyLoc2D::Type& value)
{
    static std::map<NodeBdyLoc2D::Type, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(NodeBdyLoc2D::XLO_YLO);
        INSERT_ELEMENT(NodeBdyLoc2D::XHI_YLO);
        INSERT_ELEMENT(NodeBdyLoc2D::XLO_YHI);
        INSERT_ELEMENT(NodeBdyLoc2D::XHI_YHI);
#undef INSERT_ELEMENT
    }
    
    return os << strings[value];
}


/*
 * Function to print out enum EdgeBdyLoc3D::Type value as text.
 */
inline std::ostream& operator<<(std::ostream& output, const EdgeBdyLoc3D::Type value)
{
    static std::map<EdgeBdyLoc3D::Type, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(EdgeBdyLoc3D::XLO_YLO);
        INSERT_ELEMENT(EdgeBdyLoc3D::XHI_YLO);
        INSERT_ELEMENT(EdgeBdyLoc3D::XLO_YHI);
        INSERT_ELEMENT(EdgeBdyLoc3D::XHI_YHI);
        INSERT_ELEMENT(EdgeBdyLoc3D::XLO_ZLO);
        INSERT_ELEMENT(EdgeBdyLoc3D::XHI_ZLO);
        INSERT_ELEMENT(EdgeBdyLoc3D::XLO_ZHI);
        INSERT_ELEMENT(EdgeBdyLoc3D::XHI_ZHI);
        INSERT_ELEMENT(EdgeBdyLoc3D::YLO_ZLO);
        INSERT_ELEMENT(EdgeBdyLoc3D::YHI_ZLO);
        INSERT_ELEMENT(EdgeBdyLoc3D::YLO_ZHI);
        INSERT_ELEMENT(EdgeBdyLoc3D::YHI_ZHI);
#undef INSERT_ELEMENT
    }
    
    return output << strings[value];
}


/*
 * Function to print out enum NodeBdyLoc3D::Type value as text.
 */
inline std::ostream& operator<<(std::ostream& output, const NodeBdyLoc3D::Type value)
{
    static std::map<NodeBdyLoc3D::Type, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(NodeBdyLoc3D::XLO_YLO_ZLO);
        INSERT_ELEMENT(NodeBdyLoc3D::XHI_YLO_ZLO);
        INSERT_ELEMENT(NodeBdyLoc3D::XLO_YHI_ZLO);
        INSERT_ELEMENT(NodeBdyLoc3D::XHI_YHI_ZLO);
        INSERT_ELEMENT(NodeBdyLoc3D::XLO_YLO_ZHI);
        INSERT_ELEMENT(NodeBdyLoc3D::XHI_YLO_ZHI);
        INSERT_ELEMENT(NodeBdyLoc3D::XLO_YHI_ZHI);
        INSERT_ELEMENT(NodeBdyLoc3D::XHI_YHI_ZHI);
#undef INSERT_ELEMENT
    }
    
    return output << strings[value];
}


/*
 * Function to print out enum BdryCond::Type value as text.
 */
inline std::ostream& operator<<(std::ostream& output, const BdryCond::Type value)
{
    static std::map<BdryCond::Type, std::string> strings;
    
    if (strings.size() == 0)
    {
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(BdryCond::FLOW);
        INSERT_ELEMENT(BdryCond::REFLECT);
        INSERT_ELEMENT(BdryCond::DIRICHLET);
        INSERT_ELEMENT(BdryCond::NEUMANN);
        INSERT_ELEMENT(BdryCond::XFLOW);
        INSERT_ELEMENT(BdryCond::YFLOW);
        INSERT_ELEMENT(BdryCond::ZFLOW);
        INSERT_ELEMENT(BdryCond::XREFLECT);
        INSERT_ELEMENT(BdryCond::YREFLECT);
        INSERT_ELEMENT(BdryCond::ZREFLECT);
        INSERT_ELEMENT(BdryCond::XDIRICHLET);
        INSERT_ELEMENT(BdryCond::YDIRICHLET);
        INSERT_ELEMENT(BdryCond::ZDIRICHLET);
        INSERT_ELEMENT(BdryCond::XNEUMANN);
        INSERT_ELEMENT(BdryCond::YNEUMANN);
        INSERT_ELEMENT(BdryCond::ZNEUMANN);
#undef INSERT_ELEMENT
    }
    
    return output << strings[value];
}


class EulerBoundaryConditions:
    public appu::BoundaryUtilityStrategy
{
    public:
        EulerBoundaryConditions(
            const std::string& object_name,
            const std::string& project_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& boundary_conditions_db,
            const bool& is_from_restart);
        
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
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
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
            const boost::shared_ptr<tbox::Database>& db,
            std::string& db_name,
            int bdry_location_index);
        
        /*
         * This routine is a concrete implementation of the virtual function
         * in the base class BoundaryUtilityStrategy.  It is a blank implementation
         * for the purposes of this class.
         */
        void readNeumannBoundaryDataEntry(
            const boost::shared_ptr<tbox::Database>& db,
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
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
    private:
        std::vector<double>
        readPrimitiveDataEntry(
            boost::shared_ptr<tbox::Database> db,
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
         * boost::shared_ptr to the grid geometry.
         */
        const boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Flow model.
         */
        const boost::shared_ptr<FlowModel> d_flow_model;
        
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
        
};

#endif /* EULER_BOUNDARY_CONDITIONS_HPP */
