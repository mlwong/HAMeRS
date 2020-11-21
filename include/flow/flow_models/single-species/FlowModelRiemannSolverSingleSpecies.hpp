#ifndef FLOW_MODEL_RIEMANN_SOLVER_SINGLE_SPECIES_HPP
#define FLOW_MODEL_RIEMANN_SOLVER_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModelRiemannSolver.hpp"

class FlowModelRiemannSolverSingleSpecies: public FlowModelRiemannSolver
{
    public:
        FlowModelRiemannSolverSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species):
                FlowModelRiemannSolver(
                    object_name,
                    dim,
                    grid_geometry,
                    num_species)
        {}
        
        ~FlowModelRiemannSolverSingleSpecies() {}
        
        /*
         * Compute the convective flux from conservative variables.
         */
        void
        computeConvectiveFluxFromConservativeVariables(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_minus,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type,
            const hier::Box& domain) const;
        
        /*
         * Compute the convective flux from primitive variables.
         */
        void
        computeConvectiveFluxFromPrimitiveVariables(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_minus,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type,
            const hier::Box& domain) const;
        
        /*
         * Compute the convective flux and velocity from conservative variables.
         */
        void
        computeConvectiveFluxAndVelocityFromConservativeVariables(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_minus,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type,
            const hier::Box& domain) const;
        
        /*
         * Compute the convective flux and velocity from primitive variables.
         */
        void
        computeConvectiveFluxAndVelocityFromPrimitiveVariables(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_minus,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_plus,
            const DIRECTION::TYPE& direction,
            const RIEMANN_SOLVER::TYPE& riemann_solver_type,
            const hier::Box& domain) const;
        
    private:
        /*
         * Compute the convective flux and velocity in the x-direction from conservative variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_L,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_R,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the y-direction from conservative variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_B,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_T,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the z-direction from conservative variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_B,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_F,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the x-direction from primitive variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInXDirectionFromPrimitiveVariablesHLLC(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_L,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_R,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the y-direction from primitive variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_B,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_T,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the z-direction from primitive variables with
         * HLLC Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_B,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_F,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the x-direction from conservative variables with
         * HLLC-HLL Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInXDirectionFromConservativeVariablesHLLC_HLL(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_L,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_R,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the y-direction from conservative variables with
         * HLLC-HLL Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInYDirectionFromConservativeVariablesHLLC_HLL(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_B,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_T,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the z-direction from conservative variables with
         * HLLC-HLL Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInZDirectionFromConservativeVariablesHLLC_HLL(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_B,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables_F,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the x-direction from primitive variables with
         * HLLC-HLL Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInXDirectionFromPrimitiveVariablesHLLC_HLL(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_L,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_R,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the y-direction from primitive variables with
         * HLLC-HLL Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInYDirectionFromPrimitiveVariablesHLLC_HLL(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_B,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_T,
            const hier::Box& domain,
            bool compute_velocity) const;
        
        /*
         * Compute the convective flux and velocity in the z-direction from primitive variables with
         * HLLC-HLL Riemann solver.
         */
        void
        computeConvectiveFluxAndVelocityInZDirectionFromPrimitiveVariablesHLLC_HLL(
            HAMERS_SHARED_PTR<pdat::SideData<double> > convective_flux,
            HAMERS_SHARED_PTR<pdat::SideData<double> > velocity,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_B,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables_F,
            const hier::Box& domain,
            bool compute_velocity) const;
        
};

#endif /* FLOW_MODEL_RIEMANN_SOLVER_SINGLE_SPECIES_HPP */
