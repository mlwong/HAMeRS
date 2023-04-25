#ifndef FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_SINGLE_SPECIES_HPP
#define FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModelImmersedBoundaryMethod.hpp"

class FlowModelImmersedBoundaryMethodSingleSpecies: public FlowModelImmersedBoundaryMethod
{
    public:
        FlowModelImmersedBoundaryMethodSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<ImmersedBoundaries>& immersed_boundaries,
            const HAMERS_SHARED_PTR<tbox::Database>& immersed_boundary_method_db,
            const HAMERS_SHARED_PTR<EquationOfStateMixingRules>& equation_of_state_mixing_rules);
        
        ~FlowModelImmersedBoundaryMethodSingleSpecies() {}
        
        /*
         * Set the immersed boundary method ghost cells for the cell data of conservative variables.
         */
        void setConservativeVariablesCellDataImmersedBoundaryGhosts(
            const hier::Patch& patch,
            const double data_time,
            const bool initial_time,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_var_data,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_mask,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_wall_distance,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_surface_normal,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_ip_index, // AFK
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_ip_corr, // AFK
            const hier::IntVector& offset_cons_var,
            const hier::IntVector& offset_IB,
            const hier::IntVector& ghostcell_dims_cons_var,
            const hier::IntVector& ghostcell_dims_IB,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims);
        
    private:
        /* 
         * Values of primitive variables inside the body.
         */
        Real d_rho_body;
        std::vector<Real> d_vel_body;
        Real d_p_body;
        // Real d_type;
        std::string d_type;
        
        /* 
         * Values of conservative variables inside the body.
         */
        std::vector<Real> d_mom_body;
        Real d_E_body;
        
        /*
         * Thermodynamic properties of the species.
         */
        std::vector<Real> d_thermo_properties;
        
};

#endif /* FLOW_MODEL_IMMERSED_BOUNDARY_METHOD_SINGLE_SPECIES_HPP */
