#ifndef FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_HPP
#define FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "algs/integrator/RungeKuttaLevelIntegrator.hpp"
#include "extn/visit_data_writer/ExtendedVisItDataWriter.hpp"
#include "flow/flow_models/FlowModel.hpp"
#include "util/immersed_boundaries/ImmersedBoundaries.hpp"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"

class FlowModel;

class FlowModelImmersedBoundaryMethod
{
    public:
        FlowModelImmersedBoundaryMethod(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<ImmersedBoundaries>& immersed_boundaries,
            const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db,
            const HAMERS_SHARED_PTR<EquationOfStateMixingRules>& equation_of_state_mixing_rules);
        
        virtual ~FlowModelImmersedBoundaryMethod() {}
        
        /*
         * Set the weak pointer to the flow model from the parent FlowModel class.
         */
        void setFlowModel(const HAMERS_WEAK_PTR<FlowModel>& flow_model)
        {
            d_flow_model = flow_model;
        }
        
        /*
         * Get the pointer to the immersed boundaries.
         */
        HAMERS_SHARED_PTR<ImmersedBoundaries> getImmersedBoundaries() const
        {
            return d_immersed_boundaries;
        }
        
        void putToRestart(const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db) const
        {
        }
        
        /*
         * Register the immersed boundary method variables.
         */
        void registerImmersedBoundaryMethodVariables(
            RungeKuttaLevelIntegrator* integrator,
            const hier::IntVector& num_ghosts,
            const hier::IntVector& num_ghosts_intermediate);
        
        /*
         * Register the plotting quantities.
         */
#ifdef HAVE_HDF5
        void
        registerPlotQuantities(
            const HAMERS_SHARED_PTR<ExtendedVisItDataWriter>& visit_writer,
            const HAMERS_SHARED_PTR<hier::VariableContext>& plot_context);
#endif
        
        /*
         * Set the immersed boundary method variables.
         */
        void setImmersedBoundaryMethodVariables(
            const hier::Box& domain,
            const double data_time,
            const bool initial_time,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Set the immersed boundary method ghost cells for the cell data of conservative variables.
         */
        void setConservativeVariablesCellDataImmersedBoundaryGhosts(
            const hier::Box& domain,
            const double data_time,
            const bool initial_time,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context_IB);
        
        /*
         * Set the immersed boundary method ghost cells for the cell data of conservative variables.
         */
        virtual void setConservativeVariablesCellDataImmersedBoundaryGhosts(
            const hier::Patch& patch,
            const double data_time,
            const bool initial_time,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_var_data,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_mask,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_wall_distance,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_surface_normal,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_image_points, // AFK 03/14/23
            const hier::IntVector& offset_cons_var,
            const hier::IntVector& offset_IB,
            const hier::IntVector& ghostcell_dims_cons_var,
            const hier::IntVector& ghostcell_dims_IB,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims) = 0;
        
        /*
         * Get the cell data of the immersed boundary mask in the registered patch.
         */
        HAMERS_SHARED_PTR<pdat::CellData<int> >
        getCellDataOfImmersedBoundaryMask(
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
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
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * Pointer to immersed boundaries.
         */
        HAMERS_SHARED_PTR<ImmersedBoundaries> d_immersed_boundaries;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfStateMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfStateMixingRules> d_equation_of_state_mixing_rules;
        
        /*
         * HAMERS_WEAK_PTR to FlowModel.
         */
        HAMERS_WEAK_PTR<FlowModel> d_flow_model;
        
        /*
         * HAMERS_SHARED_PTR to registered cell variables of immersed boundary methods.
         */
        static HAMERS_SHARED_PTR<pdat::CellVariable<int> > s_variable_mask;
        static HAMERS_SHARED_PTR<pdat::CellVariable<Real> > s_variable_wall_distance;
        static HAMERS_SHARED_PTR<pdat::CellVariable<Real> > s_variable_surface_normal;
        
};

#endif /* FLOW_MODEL_BASIC_UTILITIES_HPP */
