#ifndef FLOW_MODEL_BOUNDARY_UTILITIES_FIVE_EQN_ALLAIRE_HPP
#define FLOW_MODEL_BOUNDARY_UTILITIES_FIVE_EQN_ALLAIRE_HPP

#include "flow/flow_models/FlowModelBoundaryUtilities.hpp"

class FlowModelBoundaryUtilitiesFiveEqnAllaire: public FlowModelBoundaryUtilities
{
    public:
        FlowModelBoundaryUtilitiesFiveEqnAllaire(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<EquationOfStateMixingRules>& equation_of_state_mixing_rules):
                FlowModelBoundaryUtilities(
                    object_name,
                    dim,
                    num_species,
                    num_eqn,
                    equation_of_state_mixing_rules)
        {}
        
        ~FlowModelBoundaryUtilitiesFiveEqnAllaire() {}
        
        /*
         * Function to read 1d boundary data from input database.
         * Node locations that have boundary conditions identified are removed from the container.
         */
        void
        getFromInput1d(
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
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
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
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
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            std::vector<int>& face_locs,
            std::vector<int>& edge_locs,
            std::vector<int>& node_locs,
            std::vector<int>& face_conds,
            std::vector<int>& edge_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
        /*
         * Function to fill 1d node boundary values for a patch.
         * Node locations that have boundary conditions identified are removed from the container.
         */
        void
        fill1dNodeBoundaryData(
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_var_data,
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
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_var_data,
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
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_var_data,
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
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_var_data,
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
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_var_data,
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
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_node_locs,
            const std::vector<int>& bdry_node_conds,
            const std::vector<std::vector<double> >& bdry_face_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(3)));
        
    private:
        void
        read1dBdryNodes(
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            std::vector<int>& node_locs,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
        void
        read2dBdryEdges(
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            std::vector<int>& edge_locs,
            std::vector<int>& edge_conds,
            const hier::IntVector& periodic);
        
        void
        read2dBdryNodes(
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            std::vector<int>& node_locs,
            const std::vector<int>& edge_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
        void
        read3dBdryFaces(
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            std::vector<int>& face_locs,
            std::vector<int>& face_conds,
            const hier::IntVector& periodic);
        
        void
        read3dBdryEdges(
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            std::vector<int>& edge_locs,
            const std::vector<int>& face_conds,
            std::vector<int>& edge_conds,
            const hier::IntVector& periodic);
        
        void
        read3dBdryNodes(
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            std::vector<int>& node_locs,
            const std::vector<int>& face_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic);
        
};

#endif /* FLOW_MODEL_BOUNDARY_UTILITIES_FIVE_EQN_ALLAIRE_HPP */
