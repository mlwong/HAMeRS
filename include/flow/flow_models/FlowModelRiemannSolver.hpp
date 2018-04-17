#ifndef FLOW_MODEL_RIEMANN_SOLVER_HPP
#define FLOW_MODEL_RIEMANN_SOLVER_HPP

#include "HAMeRS_config.hpp"

#include "flow/flow_models/FlowModel.hpp"
#include "util/Directions.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/SideData.h"

#include "boost/weak_ptr.hpp"
#include <string>

namespace RIEMANN_SOLVER
{
    enum TYPE { HLLC,
                HLLC_HLL };
}

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
        void
        computeConvectiveFluxFromConservativeVariables(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_minus,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type) const
        {
            const hier::Box empty_box(d_dim);
            computeConvectiveFluxFromConservativeVariables(
                convective_flux,
                conservative_variables_minus,
                conservative_variables_plus,
                direction,
                riemann_solver_type,
                empty_box);
        }
        
        /*
         * Compute the convective flux from conservative variables.
         */
        virtual void
        computeConvectiveFluxFromConservativeVariables(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_minus,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the convective flux from primitive variables.
         */
        void
        computeConvectiveFluxFromPrimitiveVariables(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_minus,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type) const
        {
            const hier::Box empty_box(d_dim);
            computeConvectiveFluxFromPrimitiveVariables(
                convective_flux,
                primitive_variables_minus,
                primitive_variables_plus,
                direction,
                riemann_solver_type,
                empty_box);
        }
        
        /*
         * Compute the convective flux from primitive variables.
         */
        virtual void
        computeConvectiveFluxFromPrimitiveVariables(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_minus,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the convective flux and velocity from conservative variables.
         */
        void
        computeConvectiveFluxAndVelocityFromConservativeVariables(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            boost::shared_ptr<pdat::SideData<double> > velocity,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_minus,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type) const
        {
            const hier::Box empty_box(d_dim);
            computeConvectiveFluxAndVelocityFromConservativeVariables(
                convective_flux,
                velocity,
                conservative_variables_minus,
                conservative_variables_plus,
                direction,
                riemann_solver_type,
                empty_box);
        }
        
        /*
         * Compute the convective flux and velocity from conservative variables.
         */
        virtual void
        computeConvectiveFluxAndVelocityFromConservativeVariables(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            boost::shared_ptr<pdat::SideData<double> > velocity,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_minus,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type,
            const hier::Box& domain) const = 0;
        
        /*
         * Compute the convective flux and velocity from primitive variables.
         */
        void
        computeConvectiveFluxAndVelocityFromPrimitiveVariables(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            boost::shared_ptr<pdat::SideData<double> > velocity,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_minus,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type) const
        {
            const hier::Box empty_box(d_dim);
            computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux,
                velocity,
                primitive_variables_minus,
                primitive_variables_plus,
                direction,
                riemann_solver_type,
                empty_box);
        }
        
        /*
         * Compute the convective flux and velocity from primitive variables.
         */
        virtual void
        computeConvectiveFluxAndVelocityFromPrimitiveVariables(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            boost::shared_ptr<pdat::SideData<double> > velocity,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_minus,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type,
            const hier::Box& domain) const = 0;
        
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
