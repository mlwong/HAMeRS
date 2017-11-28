#ifndef FLOW_MODEL_RIEMANN_SOLVER_SINGLE_SPECIES_HPP
#define FLOW_MODEL_RIEMANN_SOLVER_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModelRiemannSolver.hpp"

class FlowModelRiemannSolverSingleSpecies: public FlowModelRiemannSolver
{
    public:
        FlowModelRiemannSolverSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species):
                FlowModelRiemannSolver(
                    d_object_name,
                    d_dim,
                    d_grid_geometry,
                    d_num_species)
        {}
        
        /*
         * Compute the convective flux from conservative variables.
         */
        void
        computeConvectiveFluxFromConservativeVariables(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_minus,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type,
            const hier::Box& domain) const;
        
        /*
         * Compute the convective flux from primitive variables.
         */
        void
        computeConvectiveFluxFromPrimitiveVariables(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_minus,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type,
            const hier::Box& domain) const;
        
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
            const RIEMANN_SOLVER::TYPE& riemann_solver_type,
            const hier::Box& domain) const;
        
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
            const RIEMANN_SOLVER::TYPE& riemann_solver_type,
            const hier::Box& domain) const;
        
    private:
        /*
         * Compute the convective flux and velocity in the x-direction from conservative variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC_old(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            boost::shared_ptr<pdat::SideData<double> > velocity,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_L,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_R,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the x-direction from conservative variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            boost::shared_ptr<pdat::SideData<double> > velocity,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_L,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_R,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the y-direction from conservative variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            boost::shared_ptr<pdat::SideData<double> > velocity,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_B,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_T,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the z-direction from conservative variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            boost::shared_ptr<pdat::SideData<double> > velocity,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_B,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables_F,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the x-direction from primitive variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInXDirectionFromPrimitiveVariablesHLLC(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            boost::shared_ptr<pdat::SideData<double> > velocity,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_L,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_R,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the y-direction from primitive variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            boost::shared_ptr<pdat::SideData<double> > velocity,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_B,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_T,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the z-direction from primitive variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC(
            boost::shared_ptr<pdat::SideData<double> > convective_flux,
            boost::shared_ptr<pdat::SideData<double> > velocity,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_B,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables_F,
            const hier::Box& domain,
            bool compute_velocity) const;
        
};

#endif /* FLOW_MODEL_RIEMANN_SOLVER_SINGLE_SPECIES_HPP */
