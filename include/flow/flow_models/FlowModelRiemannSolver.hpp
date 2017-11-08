#ifndef FLOW_MODEL_RIEMANN_SOLVER_HPP
#define FLOW_MODEL_RIEMANN_SOLVER_HPP

#include "HAMeRS_config.hpp"

#include "flow/flow_models/FlowModel.hpp"
#include "flow/flow_models/RiemannSolverTypes.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/SideData.h"

#include "boost/weak_ptr.hpp"
#include <string>

class FlowModel;

class FlowModelRiemannSolver
{
    public:
        FlowModelRiemannSolver(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_species(num_species)
        {}
        
        /*
         * Set the weak pointer to the flow model from the parent FlowModel class.
         */
        void setFlowModel(const boost::weak_ptr<FlowModel>& flow_model)
        {
            d_flow_model = flow_model;
        }
        
        /*
         * Compute the convective flux from conservative variables.
         */
        virtual void
        computeConvectiveFluxFromConservativeVariables(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
            const hier::Box& domain,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type) = 0;
        
        /*
         * Compute the convective flux from primitive variables.
         */
        virtual void
        computeConvectiveFluxFromPrimitiveVariables(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables,
            const hier::Box& domain,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type) = 0;
    
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
         * boost::shared_ptr to the grid geometry.
         */
        const boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * boost::weak_ptr to FlowModel.
         */
        boost::weak_ptr<FlowModel> d_flow_model;
        
};

#endif /* FLOW_MODEL_RIEMANN_SOLVER_HPP */
