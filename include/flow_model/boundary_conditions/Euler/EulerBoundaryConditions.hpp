#ifndef EULER_BOUNDARY_CONDITIONS_HPP
#define EULER_BOUNDARY_CONDITIONS_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/appu/CartesianBoundaryDefines.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Dimension.h"

#include "equation_of_state/EquationOfStateIdealGas.hpp"
#include "flow_model/FlowModels.hpp"

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
            const FLOW_MODEL& flow_model,
            const int& num_species,
            const boost::shared_ptr<EquationOfState>& equation_of_state,
            const boost::shared_ptr<tbox::Database>& boundary_conditions_db,
            const bool& is_from_restart);
        
        /*
         * Set the number of ghost cells needed.
         */
        void
        setNumberOfGhostCells(const hier::IntVector& num_ghosts)
        {
            d_num_ghosts = num_ghosts;
            
            d_num_ghosts_set = true;
        }
        
        /*
         * Set the cell variables if single-species flow model is chosen.
         */
        void
        setVariablesForSingleSpecies(
            const boost::shared_ptr<pdat::CellVariable<double> >& density,
            const boost::shared_ptr<pdat::CellVariable<double> >& momentum,
            const boost::shared_ptr<pdat::CellVariable<double> >& total_energy)
        {
            if (d_flow_model != SINGLE_SPECIES)
            {            
                TBOX_ERROR(d_object_name
                           << ": "
                           << "setVariablesForSingleSpecies() shouldn't be used."
                           << std::endl);
            }
            
            d_density = density;
            d_momentum = momentum;
            d_total_energy = total_energy;
            
            d_variables_set = true;
        }
        
        /*
         * Set the cell variables if conservative four-equation multi-species flow model
         * is chosen.
         */
        void
        setVariablesForFourEqnConservative(
            const boost::shared_ptr<pdat::CellVariable<double> > partial_density,
            const boost::shared_ptr<pdat::CellVariable<double> > momentum,
            const boost::shared_ptr<pdat::CellVariable<double> > total_energy)
        {
            if (d_flow_model != FOUR_EQN_CONSERVATIVE)
            {            
                TBOX_ERROR(d_object_name
                           << ": "
                           << "setVariablesForFourEqnConservative() shouldn't be used."
                           << std::endl);
            }
            
            d_partial_density = partial_density;
            d_momentum = momentum;
            d_total_energy = total_energy;
            
            d_variables_set = true;
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
            if (d_flow_model != FOUR_EQN_SHYUE)
            {            
                TBOX_ERROR(d_object_name
                           << ": "
                           << "setVariablesForFourEqnShyue() shouldn't be used."
                           << std::endl);
            }
            
            d_density = density;
            d_momentum = momentum;
            d_total_energy = total_energy;
            d_mass_fraction = mass_fraction;
            
            d_variables_set = true;
        }
        
        /*
         * Set the cell variables if five-equation multi-species flow model
         * by Allaire is chosen.
         */
        void
        setVariablesForFiveEqnAllaire(
            const boost::shared_ptr<pdat::CellVariable<double> > partial_density,
            const boost::shared_ptr<pdat::CellVariable<double> > momentum,
            const boost::shared_ptr<pdat::CellVariable<double> > total_energy,
            const boost::shared_ptr<pdat::CellVariable<double> > volume_fraction)
        {
            if (d_flow_model != FIVE_EQN_ALLAIRE)
            {            
                TBOX_ERROR(d_object_name
                           << ": "
                           << "setVariablesForFiveEqnAllaire() shouldn't be used."
                           << std::endl);
            }
            
            d_partial_density = partial_density;
            d_momentum = momentum;
            d_total_energy = total_energy;
            d_volume_fraction = volume_fraction;
            
            d_variables_set = true;
        }
        
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
        void
        readStateDataEntryForSingleSpecies(
            boost::shared_ptr<tbox::Database> db,
            const std::string& db_name,
            int array_indx,
            std::vector<double>& density,
            std::vector<double>& momentum,
            std::vector<double>& total_energy);
        
        void
        readStateDataEntryForFourEqnConservative(
            boost::shared_ptr<tbox::Database> db,
            const std::string& db_name,
            int array_indx,
            std::vector<double>& partial_density,
            std::vector<double>& momentum,
            std::vector<double>& total_energy);
        
        void
        readStateDataEntryForFourEqnShyue(
            boost::shared_ptr<tbox::Database> db,
            const std::string& db_name,
            int array_indx,
            std::vector<double>& density,
            std::vector<double>& momentum,
            std::vector<double>& total_energy,
            std::vector<double>& mass_fraction);
        
        void
        readStateDataEntryForFiveEqnAllaire(
            boost::shared_ptr<tbox::Database> db,
            const std::string& db_name,
            int array_indx,
            std::vector<double>& partial_density,
            std::vector<double>& momentum,
            std::vector<double>& total_energy,
            std::vector<double>& volume_fraction);
        
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
         * Number of ghost cells for time-independent variables.
         */
        hier::IntVector d_num_ghosts;
        
        /*
         * Flow model.
         */
        const FLOW_MODEL d_flow_model;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * boost::shared_ptr to EquationOfState.
         */
        const boost::shared_ptr<EquationOfState> d_equation_of_state;
        
        /*
         * boost::shared_ptr to solution state.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_partial_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_momentum;
        boost::shared_ptr<pdat::CellVariable<double> > d_total_energy;
        boost::shared_ptr<pdat::CellVariable<double> > d_mass_fraction;
        boost::shared_ptr<pdat::CellVariable<double> > d_volume_fraction;
        
        /*
         * Boolean to determine whether proper variables are initialized.
         */
        bool d_variables_set;
        
        /*
         * Boolean to determine whether the number of ghost cells is initialized.
         */
        bool d_num_ghosts_set;
        
        /*
         * Boundary condition cases and boundary values.
         * Options are: FLOW, REFLECT, DIRICHLET
         * and variants for nodes and edges.
         *
         * Input file values are read into these arrays.
         */
        std::vector<int> d_master_bdry_node_conds;
        std::vector<int> d_master_bdry_edge_conds; // Used in 2D and 3D only.
        std::vector<int> d_master_bdry_face_conds; // Used in 3D only.
        
        /*
         * Boundary condition cases for scalar and vector (i.e., depth > 1)
         * variables.  These are post-processed input values and are passed
         * to the boundary routines.
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
        std::vector<double> d_bdry_node_density;                // Used in 1D only.
        std::vector<double> d_bdry_node_partial_density;        // Used in 1D only.
        std::vector<double> d_bdry_node_momentum;               // Used in 1D only.
        std::vector<double> d_bdry_node_total_energy;           // Used in 1D only.
        std::vector<double> d_bdry_node_mass_fraction;          // Used in 1D only.
        std::vector<double> d_bdry_node_volume_fraction;        // Used in 1D only.
        
        std::vector<double> d_bdry_edge_density;                // Used in 2D only.
        std::vector<double> d_bdry_edge_partial_density;        // Used in 2D only.
        std::vector<double> d_bdry_edge_momentum;               // Used in 2D only.
        std::vector<double> d_bdry_edge_total_energy;           // Used in 2D only.
        std::vector<double> d_bdry_edge_mass_fraction;          // Used in 2D only.
        std::vector<double> d_bdry_edge_volume_fraction;        // Used in 2D only.
        
        std::vector<double> d_bdry_face_density;                // Used in 3D only.
        std::vector<double> d_bdry_face_partial_density;        // Used in 3D only.
        std::vector<double> d_bdry_face_momentum;               // Used in 3D only.
        std::vector<double> d_bdry_face_total_energy;           // Used in 3D only.
        std::vector<double> d_bdry_face_mass_fraction;          // Used in 3D only.
        std::vector<double> d_bdry_face_volume_fraction;        // Used in 3D only.
    
};

#endif /* EULER_BOUNDARY_CONDITIONS_HPP */
