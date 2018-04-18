#ifndef FLOW_MODEL_BOUNDARY_UTILITIES_HPP
#define FLOW_MODEL_BOUNDARY_UTILITIES_HPP

#include "HAMeRS_config.hpp"

#include "util/mixing_rules/equations_of_state/EquationOfStateMixingRulesManager.hpp"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Database.h"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

namespace BDRY_COND
{
    namespace FLOW_MODEL
    {
        enum TYPE
        {
            ADIABATIC_NO_SLIP   = 10,
            ISOTHERMAL_NO_SLIP  = 11,
            XADIABATIC_NO_SLIP  = 100,
            YADIABATIC_NO_SLIP  = 101,
            ZADIABATIC_NO_SLIP  = 102,
            XISOTHERMAL_NO_SLIP = 110,
            YISOTHERMAL_NO_SLIP = 111,
            ZISOTHERMAL_NO_SLIP = 112
        };
    }
}

class FlowModelBoundaryUtilities
{
    public:
        FlowModelBoundaryUtilities(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const int& num_species,
            const int& num_eqn,
            const boost::shared_ptr<EquationOfStateMixingRules>& equation_of_state_mixing_rules):
                d_object_name(object_name),
                d_dim(dim),
                d_num_species(num_species),
                d_num_eqn(num_eqn),
                d_equation_of_state_mixing_rules(equation_of_state_mixing_rules)
        {}
        
        virtual ~FlowModelBoundaryUtilities() {}
        
        /*
         * Virtual function to read 1d boundary data from input database.
         * Node locations that have boundary conditions identified are removed from the container.
         */
        virtual void
        getFromInput1d(
            const boost::shared_ptr<tbox::Database>& input_db,
            std::vector<int>& node_locs,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic) = 0;
        
        /*
         * Virtual function to read 2d boundary data from input database.
         * Node and edge locations that have boundary conditions identified are removed from the
         * containers.
         */
        virtual void
        getFromInput2d(
            const boost::shared_ptr<tbox::Database>& input_db,
            std::vector<int>& edge_locs,
            std::vector<int>& node_locs,
            std::vector<int>& edge_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic) = 0;
        
        /*
         * Virtual function to read 3d boundary data from input database.
         * Node, edge and face locations that have boundary conditions identified are removed from
         * the containers.
         */
        virtual void
        getFromInput3d(
            const boost::shared_ptr<tbox::Database>& input_db,
            std::vector<int>& face_locs,
            std::vector<int>& edge_locs,
            std::vector<int>& node_locs,
            std::vector<int>& face_conds,
            std::vector<int>& edge_conds,
            std::vector<int>& node_conds,
            const hier::IntVector& periodic) = 0;
        
        /*
         * Virtual unction to fill 1d node boundary values for a patch.
         * Node locations that have boundary conditions identified are removed from the container.
         */
        virtual void
        fill1dNodeBoundaryData(
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_node_locs,
            const std::vector<int>& bdry_node_conds,
            const std::vector<std::vector<double> >& bdry_node_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(1))) = 0;
        
        /*
         * Virtual function to fill 2d edge boundary values for a patch.
         * Edge locations that have boundary conditions identified are removed from the container.
         */
        virtual void
        fill2dEdgeBoundaryData(
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_edge_locs,
            const std::vector<int>& bdry_edge_conds,
            const std::vector<std::vector<double> >& bdry_edge_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(2))) = 0;
        
        /*
         * Virtual function to fill 2d node boundary values for a patch.
         * Node locations that have boundary conditions identified are removed from the container.
         */
        virtual void
        fill2dNodeBoundaryData(
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_node_locs,
            const std::vector<int>& bdry_node_conds,
            const std::vector<std::vector<double> >& bdry_edge_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(2))) = 0;
        
        /*
         * Virtual function to fill 3d face boundary values for a patch.
         * Face locations that have boundary conditions identified are removed from the container.
         */
        virtual void
        fill3dFaceBoundaryData(
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_face_locs,
            const std::vector<int>& bdry_face_conds,
            const std::vector<std::vector<double> >& bdry_face_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(3))) = 0;
        
        /*
         * Virtual function to fill 3d edge boundary values for a patch.
         * Edge locations that have boundary conditions identified are removed from the container.
         */
        virtual void
        fill3dEdgeBoundaryData(
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_edge_locs,
            const std::vector<int>& bdry_edge_conds,
            const std::vector<std::vector<double> >& bdry_face_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(3))) = 0;
        
        /*
         * Virtual function to fill 3d node boundary values for a patch.
         * Node locations that have boundary conditions identified are removed from the container.
         */
        virtual void
        fill3dNodeBoundaryData(
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
            const hier::Patch& patch,
            std::vector<int>& bdry_node_locs,
            const std::vector<int>& bdry_node_conds,
            const std::vector<std::vector<double> >& bdry_face_values,
            const hier::IntVector& ghost_width_to_fill = -hier::IntVector::getOne(tbox::Dimension(3))) = 0;
        
protected:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * boost::shared_ptr to EquationOfStateMixingRules.
         */
        boost::shared_ptr<EquationOfStateMixingRules> d_equation_of_state_mixing_rules;
        
};

#endif /* FLOW_MODEL_BOUNDARY_UTILITIES_HPP */
