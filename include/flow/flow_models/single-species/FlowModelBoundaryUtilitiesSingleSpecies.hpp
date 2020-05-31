#ifndef FLOW_MODEL_BOUNDARY_UTILITIES_SINGLE_SPECIES_HPP
#define FLOW_MODEL_BOUNDARY_UTILITIES_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModelBoundaryUtilities.hpp"

class FlowModelBoundaryUtilitiesSingleSpecies: public FlowModelBoundaryUtilities
{
    public:
        FlowModelBoundaryUtilitiesSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const int& num_eqn,
            const boost::shared_ptr<EquationOfStateMixingRules>& equation_of_state_mixing_rules);
        
        ~FlowModelBoundaryUtilitiesSingleSpecies() {}
        
        /*
         * Function to read 1d boundary data from input database.
         * Node locations that have boundary conditions identified are removed from the container.
         */
        void
        getFromInput1d(
            const boost::shared_ptr<tbox::Database>& input_db,
            std::vector<int>& node_locs,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
        /*
         * Function to read 2d boundary data from input database.
         * Node and edge locations that have boundary conditions identified are removed from the
         * containers.
         */
        void
        getFromInput2d(
            const boost::shared_ptr<tbox::Database>& input_db,
            std::vector<int>& edge_locs,
            std::vector<int>& node_locs,
            std::vector<int>& edge_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
        /*
         * Function to read 3d boundary data from input database.
         * Node, edge and face locations that have boundary conditions identified are removed from
         * the containers.
         */
        void
        getFromInput3d(
            const boost::shared_ptr<tbox::Database>& input_db,
            std::vector<int>& face_locs,
            std::vector<int>& edge_locs,
            std::vector<int>& node_locs,
            std::vector<int>& face_conds,
            std::vector<int>& edge_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
        /*
         * Function that returns the integer edge boundary location
         * corresponding to the given node location and node boundary
         * condition.
         *
         * If the node boundary condition type or node location are unknown,
         * or the boundary condition type is inconsistent with the node location,
         * an error code (-1) is returned.
         */
        int
        getEdgeLocationForNodeBdry(
            int node_loc,
            int node_btype);
        
        /*
         * Function that returns the integer face boundary location
         * corresponding to the given edge location and edge boundary
         * condition.
         *
         * If the edge boundary condition type or edge location are unknown,
         * or the boundary condition type is inconsistent with the edge location,
         * an error code (-1) is returned.
         */
        int
        getFaceLocationForEdgeBdry(
            int edge_loc,
            int edge_btype);
        
        /*
         * Function that returns the integer face boundary location
         * corresponding to the given node location and node boundary
         * condition.
         *
         * If the node boundary condition type or node location are unknown,
         * or the boundary condition type is inconsistent with the node location,
         * an error code (-1) is returned.
         */
        int
        getFaceLocationForNodeBdry(
            int node_loc,
            int node_btype);
        
        /*
         * Function to fill 1d node boundary values for a patch.
         * Node locations that have boundary conditions identified are removed from the container.
         */
        void
        fill1dNodeBoundaryData(
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_node_locs,
            const std::vector<int>& bdry_node_conds,
            const std::vector<std::vector<double> >& bdry_node_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(1)));
        
        /*
         * Function to fill 2d edge boundary values for a patch.
         * Edge locations that have boundary conditions identified are removed from the container.
         */
        void
        fill2dEdgeBoundaryData(
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_edge_locs,
            const std::vector<int>& bdry_edge_conds,
            const std::vector<std::vector<double> >& bdry_edge_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(2)));
        
        /*
         * Function to fill 2d node boundary values for a patch.
         * Node locations that have boundary conditions identified are removed from the container.
         */
        void
        fill2dNodeBoundaryData(
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_node_locs,
            const std::vector<int>& bdry_node_conds,
            const std::vector<std::vector<double> >& bdry_edge_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(2)));
        
        /*
         * Function to fill 3d face boundary values for a patch.
         * Face locations that have boundary conditions identified are removed from the container.
         */
        void
        fill3dFaceBoundaryData(
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_face_locs,
            const std::vector<int>& bdry_face_conds,
            const std::vector<std::vector<double> >& bdry_face_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(3)));
        
        /*
         * Function to fill 3d edge boundary values for a patch.
         * Edge locations that have boundary conditions identified are removed from the container.
         */
        void
        fill3dEdgeBoundaryData(
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_edge_locs,
            const std::vector<int>& bdry_edge_conds,
            const std::vector<std::vector<double> >& bdry_face_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(3)));
        
        /*
         * Function to fill 3d node boundary values for a patch.
         * Node locations that have boundary conditions identified are removed from the container.
         */
        void
        fill3dNodeBoundaryData(
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_node_locs,
            const std::vector<int>& bdry_node_conds,
            const std::vector<std::vector<double> >& bdry_face_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(3)));
        
    private:
        void
        read1dBdryNodes(
            const boost::shared_ptr<tbox::Database>& input_db,
            std::vector<int>& node_locs,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
        void
        read2dBdryEdges(
            const boost::shared_ptr<tbox::Database>& input_db,
            std::vector<int>& edge_locs,
            std::vector<int>& edge_conds,
            const hier::IntVector& periodic);
        
        void
        read2dBdryNodes(
            const boost::shared_ptr<tbox::Database>& input_db,
            std::vector<int>& node_locs,
            const std::vector<int>& edge_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
        void
        read3dBdryFaces(
            const boost::shared_ptr<tbox::Database>& input_db,
            std::vector<int>& face_locs,
            std::vector<int>& face_conds,
            const hier::IntVector& periodic);
        
        void
        read3dBdryEdges(
            const boost::shared_ptr<tbox::Database>& input_db,
            std::vector<int>& edge_locs,
            const std::vector<int>& face_conds,
            std::vector<int>& edge_conds,
            const hier::IntVector& periodic);
        
        void
        read3dBdryNodes(
            const boost::shared_ptr<tbox::Database>& input_db,
            std::vector<int>& node_locs,
            const std::vector<int>& face_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
        void
        readAdiabaticNoSlip(
            const boost::shared_ptr<tbox::Database>& db,
            std::string& db_name,
            int bdry_location_index);
        
        void
        readIsothermalNoSlip(
            const boost::shared_ptr<tbox::Database>& db,
            std::string& db_name,
            int bdry_location_index);
        
        /*
         * Thermodynamic properties of the species.
         */
        std::vector<double> d_thermo_properties;
        
        /*
         * Vectors of node (1D), edge (2D) or face (3D) boundary values for ADIABATIC_NO_SLIP case.
         */
        std::vector<double> d_bdry_node_adiabatic_no_slip_vel;
        std::vector<double> d_bdry_edge_adiabatic_no_slip_vel;
        std::vector<double> d_bdry_face_adiabatic_no_slip_vel;
        
        /*
         * Vectors of node (1D), edge (2D) or face (3D) boundary values for ISOTHERMAL_NO_SLIP case.
         */
        std::vector<double> d_bdry_node_isothermal_no_slip_T;
        std::vector<double> d_bdry_edge_isothermal_no_slip_T;
        std::vector<double> d_bdry_face_isothermal_no_slip_T;
        std::vector<double> d_bdry_node_isothermal_no_slip_vel;
        std::vector<double> d_bdry_edge_isothermal_no_slip_vel;
        std::vector<double> d_bdry_face_isothermal_no_slip_vel;
        
};

#endif /* FLOW_MODEL_BOUNDARY_UTILITIES_SINGLE_SPECIES_HPP */
