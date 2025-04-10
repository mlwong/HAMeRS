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

namespace VELOCITY_IBC
{
    enum TYPE { NONE,
                SLIP,
                NO_SLIP };
}

namespace TEMPERATURE_IBC
{
    enum TYPE { NONE,
                ADIABATIC,
                ISOTHERMAL };
}

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
         * Get the additional number of ghost cells needed by the immersed boundary method.
         */
        hier::IntVector
        getImmersedBoundaryMethodAdditionalNumberOfGhostCells() const
        {
            return d_num_IBM_ghosts;
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
        
        /*
         * Compute the data on the surface triangulation.
         */
        void computeSurfaceTriangulationData(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Output the surface triangulation with surface data.
         */
        void writeSurfaceTriangulationWithData(const std::string& file_name) const;
        
    protected:
        /*
         * Dot product of two 2D vectors.
         */
        static inline __attribute__((always_inline)) Real dotProduct2D(
            const Real& a0,
            const Real& a1,
            const Real& b0,
            const Real& b1)
        {
            return a0*b0 + a1*b1;
        }
        
        /*
         * Dot product of two 3D vectors.
         */
        static inline __attribute__((always_inline)) Real dotProduct3D(
            const Real& a0,
            const Real& a1,
            const Real& a2,
            const Real& b0,
            const Real& b1,
            const Real& b2)
        {
            return a0*b0 + a1*b1 + a2*b2;
        }
        
        /*
         * Dirichlet boundary condition with second order of accuracy.
         */
        static inline __attribute__((always_inline)) Real getGhostValueDirichletBC(
            const Real& u_body,
            const Real& u_ip,
            const Real& d_ip,
            const Real& d_gc)
        {
            const Real u_gc = u_ip - ((d_ip + d_gc)/d_ip)*(u_ip - u_body);
            return u_gc;
        }
        
        /*
         * Neumann boundary condition with second order of accuracy (zero gradient).
         */
        static inline __attribute__((always_inline)) Real getGhostValueNeumannBC(
            const Real& u_ip1,
            const Real& u_ip2,
            const Real& d_ip1,
            const Real& d_ip2,
            const Real& d_gc)
        {
            const Real u_gc = (u_ip1*(d_ip2*d_ip2 - d_gc*d_gc) - u_ip2*(d_ip1*d_ip1 - d_gc*d_gc))/
                (d_ip2*d_ip2 - d_ip1*d_ip1);
            return u_gc;
        }
        
        /*
         * Neumann boundary condition with second order of accuracy (zero gradient).
         */
        static inline __attribute__((always_inline)) Real getGhostValueNeumannBC(
            const Real& dudn_body,
            const Real& u_ip1,
            const Real& u_ip2,
            const Real& d_ip1,
            const Real& d_ip2,
            const Real& d_gc)
        {
            const Real u_gc = (u_ip1*(d_ip2*d_ip2 - d_gc*d_gc) - u_ip2*(d_ip1*d_ip1 - d_gc*d_gc))/
                (d_ip2*d_ip2 - d_ip1*d_ip1) -
                (d_ip1*d_ip2 + d_gc*d_gc + d_gc*d_ip1 + d_gc*d_ip2)/(d_ip1 + d_ip2)*dudn_body;
            return u_gc;
        }
        
        /*
         * Get the indices for the 2D bilinear interpolation.
         */
         static inline __attribute__((always_inline)) void getBilinearInterpolationIndices2D(
            int& idx_BL,
            int& idx_BR,
            int& idx_TL,
            int& idx_TR,
            int& idx_IB_BL,
            int& idx_IB_BR,
            int& idx_IB_TL,
            int& idx_IB_TR,
            Real& x_ip_BL,
            Real& y_ip_BL,
            const Real& x_ip,
            const Real& y_ip,
            const Real& patch_xlo_0,
            const Real& patch_xlo_1,
            const int& offset_0,
            const int& offset_1,
            const int& ghostcell_dim_0,
            const int& offset_0_IB,
            const int& offset_1_IB,
            const int& ghostcell_dim_0_IB,
            const Real& dx,
            const Real& dx_inv)
        {
            constexpr Real half = Real(1)/Real(2);
            
            const int ip_i = int(floor((x_ip - patch_xlo_0 - half * dx)*dx_inv));
            const int ip_j = int(floor((y_ip - patch_xlo_1 - half * dx)*dx_inv));

            idx_IB_BL  = (ip_i     + offset_0_IB) + (ip_j     + offset_1_IB) * ghostcell_dim_0_IB;
            idx_IB_BR  = (ip_i + 1 + offset_0_IB) + (ip_j     + offset_1_IB) * ghostcell_dim_0_IB;
            idx_IB_TL  = (ip_i     + offset_0_IB) + (ip_j + 1 + offset_1_IB) * ghostcell_dim_0_IB;
            idx_IB_TR  = (ip_i + 1 + offset_0_IB) + (ip_j + 1 + offset_1_IB) * ghostcell_dim_0_IB;
            
            idx_BL  = (ip_i     + offset_0) + (ip_j     + offset_1) * ghostcell_dim_0;
            idx_BR  = (ip_i + 1 + offset_0) + (ip_j     + offset_1) * ghostcell_dim_0;
            idx_TL  = (ip_i     + offset_0) + (ip_j + 1 + offset_1) * ghostcell_dim_0;
            idx_TR  = (ip_i + 1 + offset_0) + (ip_j + 1 + offset_1) * ghostcell_dim_0;
            
            x_ip_BL = patch_xlo_0 + (Real(ip_i) + half)*dx;
            y_ip_BL = patch_xlo_1 + (Real(ip_j) + half)*dx;
        }
        
        /*
         * 2D bilinear interpolation.
         */
        static inline __attribute__((always_inline)) Real bilinearInterpolate2D(
            const Real& u_BL,
            const Real& u_BR,
            const Real& u_TL,
            const Real& u_TR,
            const Real& x_ip,
            const Real& y_ip,
            const Real& x_ip_BL,
            const Real& y_ip_BL,
            const Real& dx_inv)
        {
            constexpr Real one = Real(1);
            
            const Real ip_ratio_0 = (x_ip - x_ip_BL)*dx_inv;
            const Real ip_ratio_1 = (y_ip - y_ip_BL)*dx_inv;
            
            const Real u_f1 = (one - ip_ratio_0)*u_BL + ip_ratio_0*u_BR;
            const Real u_f2 = (one - ip_ratio_0)*u_TL + ip_ratio_0*u_TR;
            const Real u_ip = (one - ip_ratio_1)*u_f1 + ip_ratio_1*u_f2;
            
            return u_ip;
        }
        
        /*
         * Get the indices for the 3D trilinear interpolation.
         */
         static inline __attribute__((always_inline)) void getTrilinearInterpolationIndices3D(
            int& idx_cons_var_LBK,
            int& idx_cons_var_RBK,
            int& idx_cons_var_LTK,
            int& idx_cons_var_RTK,
            int& idx_cons_var_LBF,
            int& idx_cons_var_RBF,
            int& idx_cons_var_LTF,
            int& idx_cons_var_RTF,
            int& idx_IB_LBK,
            int& idx_IB_RBK,
            int& idx_IB_LTK,
            int& idx_IB_RTK,
            int& idx_IB_LBF,
            int& idx_IB_RBF,
            int& idx_IB_LTF,
            int& idx_IB_RTF,
            Real& x_ip_LBK,
            Real& y_ip_LBK,
            Real& z_ip_LBK,
            const Real& x_ip,
            const Real& y_ip,
            const Real& z_ip,
            const Real& patch_xlo_0,
            const Real& patch_xlo_1,
            const Real& patch_xlo_2,
            const int& offset_0,
            const int& offset_1,
            const int& offset_2,
            const int& ghostcell_dim_0,
            const int& ghostcell_dim_1,
            const int& offset_0_IB,
            const int& offset_1_IB,
            const int& offset_2_IB,
            const int& ghostcell_dim_0_IB,
            const int& ghostcell_dim_1_IB,
            const Real& dx,
            const Real& dx_inv)
        {
            constexpr Real half = Real(1)/Real(2);
            
            const int ip_i = int(floor((x_ip - patch_xlo_0 - half * dx)*dx_inv));
            const int ip_j = int(floor((y_ip - patch_xlo_1 - half * dx)*dx_inv));
            const int ip_k = int(floor((z_ip - patch_xlo_2 - half * dx)*dx_inv));

            idx_IB_LBK = (ip_i     + offset_0_IB) + (ip_j     + offset_1_IB) * ghostcell_dim_0_IB + (ip_k     + offset_2_IB) * ghostcell_dim_0_IB * ghostcell_dim_1_IB;
            idx_IB_RBK = (ip_i + 1 + offset_0_IB) + (ip_j     + offset_1_IB) * ghostcell_dim_0_IB + (ip_k     + offset_2_IB) * ghostcell_dim_0_IB * ghostcell_dim_1_IB;
            idx_IB_LTK = (ip_i     + offset_0_IB) + (ip_j + 1 + offset_1_IB) * ghostcell_dim_0_IB + (ip_k     + offset_2_IB) * ghostcell_dim_0_IB * ghostcell_dim_1_IB;
            idx_IB_RTK = (ip_i + 1 + offset_0_IB) + (ip_j + 1 + offset_1_IB) * ghostcell_dim_0_IB + (ip_k     + offset_2_IB) * ghostcell_dim_0_IB * ghostcell_dim_1_IB;
            idx_IB_LBF = (ip_i     + offset_0_IB) + (ip_j     + offset_1_IB) * ghostcell_dim_0_IB + (ip_k + 1 + offset_2_IB) * ghostcell_dim_0_IB * ghostcell_dim_1_IB;
            idx_IB_RBF = (ip_i + 1 + offset_0_IB) + (ip_j     + offset_1_IB) * ghostcell_dim_0_IB + (ip_k + 1 + offset_2_IB) * ghostcell_dim_0_IB * ghostcell_dim_1_IB;
            idx_IB_LTF = (ip_i     + offset_0_IB) + (ip_j + 1 + offset_1_IB) * ghostcell_dim_0_IB + (ip_k + 1 + offset_2_IB) * ghostcell_dim_0_IB * ghostcell_dim_1_IB;
            idx_IB_RTF = (ip_i + 1 + offset_0_IB) + (ip_j + 1 + offset_1_IB) * ghostcell_dim_0_IB + (ip_k + 1 + offset_2_IB) * ghostcell_dim_0_IB * ghostcell_dim_1_IB;
            
            idx_cons_var_LBK = (ip_i     + offset_0) + (ip_j     + offset_1) * ghostcell_dim_0 + (ip_k     + offset_2) * ghostcell_dim_0 * ghostcell_dim_1;
            idx_cons_var_RBK = (ip_i + 1 + offset_0) + (ip_j     + offset_1) * ghostcell_dim_0 + (ip_k     + offset_2) * ghostcell_dim_0 * ghostcell_dim_1;
            idx_cons_var_LTK = (ip_i     + offset_0) + (ip_j + 1 + offset_1) * ghostcell_dim_0 + (ip_k     + offset_2) * ghostcell_dim_0 * ghostcell_dim_1;
            idx_cons_var_RTK = (ip_i + 1 + offset_0) + (ip_j + 1 + offset_1) * ghostcell_dim_0 + (ip_k     + offset_2) * ghostcell_dim_0 * ghostcell_dim_1;
            idx_cons_var_LBF = (ip_i     + offset_0) + (ip_j     + offset_1) * ghostcell_dim_0 + (ip_k + 1 + offset_2) * ghostcell_dim_0 * ghostcell_dim_1;
            idx_cons_var_RBF = (ip_i + 1 + offset_0) + (ip_j     + offset_1) * ghostcell_dim_0 + (ip_k + 1 + offset_2) * ghostcell_dim_0 * ghostcell_dim_1;
            idx_cons_var_LTF = (ip_i     + offset_0) + (ip_j + 1 + offset_1) * ghostcell_dim_0 + (ip_k + 1 + offset_2) * ghostcell_dim_0 * ghostcell_dim_1;
            idx_cons_var_RTF = (ip_i + 1 + offset_0) + (ip_j + 1 + offset_1) * ghostcell_dim_0 + (ip_k + 1 + offset_2) * ghostcell_dim_0 * ghostcell_dim_1;
            
            x_ip_LBK = patch_xlo_0 + (Real(ip_i) + half)*dx;
            y_ip_LBK = patch_xlo_1 + (Real(ip_j) + half)*dx;
            z_ip_LBK = patch_xlo_2 + (Real(ip_k) + half)*dx;
        }
        
        /*
         * 3D trilinear interpolation.
         */
        static inline __attribute__((always_inline)) Real trilinearInterpolate3D(
            const Real& u_ip_LBK,
            const Real& u_ip_RBK,
            const Real& u_ip_LTK,
            const Real& u_ip_RTK,
            const Real& u_ip_LBF,
            const Real& u_ip_RBF,
            const Real& u_ip_LTF,
            const Real& u_ip_RTF,
            const Real& x_ip,
            const Real& y_ip,
            const Real& z_ip,
            const Real& x_ip_LBK,
            const Real& y_ip_LBK,
            const Real& z_ip_LBK,
            const Real& dx_inv)
        {
            constexpr Real one = Real(1);
            
            const Real ip_ratio_0 = (x_ip - x_ip_LBK)*dx_inv;
            const Real ip_ratio_1 = (y_ip - y_ip_LBK)*dx_inv;
            const Real ip_ratio_2 = (z_ip - z_ip_LBK)*dx_inv;
            
            const Real u_ip_BK = (one - ip_ratio_0)*u_ip_LBK + ip_ratio_0*u_ip_RBK;
            const Real u_ip_TK = (one - ip_ratio_0)*u_ip_LTK + ip_ratio_0*u_ip_RTK;
            const Real u_ip_BF = (one - ip_ratio_0)*u_ip_LBF + ip_ratio_0*u_ip_RBF;
            const Real u_ip_TF = (one - ip_ratio_0)*u_ip_LTF + ip_ratio_0*u_ip_RTF;
            
            const Real u_ip_F = (one - ip_ratio_1)*u_ip_BF + ip_ratio_1*u_ip_TF;
            const Real u_ip_K = (one - ip_ratio_1)*u_ip_BK + ip_ratio_1*u_ip_TK;
            
            const Real u_ip = (one - ip_ratio_2) * u_ip_K  + ip_ratio_2 * u_ip_F;
            
            return u_ip;
        }
        
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
         * Number of ghost cells needed by the immersed boundary method.
         */
        hier::IntVector d_num_IBM_ghosts;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * Type of velocity immersed boundary condition.
         */
        VELOCITY_IBC::TYPE d_bc_type_velocity;
         
        /*
         * Type of temperature immersed boundary condition.
         */
        TEMPERATURE_IBC::TYPE d_bc_type_temperature;
        
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
        
        /*
         * Data for the surface triangulation if needed.
         */
        
        std::vector<double> d_surface_triangulation_dx_grid;
        
        std::vector<std::array<double, 3> > d_surface_triangulation_coor_ip_1;
        std::vector<std::array<double, 3> > d_surface_triangulation_coor_ip_2;
        
};

#endif /* FLOW_MODEL_BASIC_UTILITIES_HPP */
