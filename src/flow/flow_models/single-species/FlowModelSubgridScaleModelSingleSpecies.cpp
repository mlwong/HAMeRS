#include "flow/flow_models/single-species/FlowModelSubgridScaleModelSingleSpecies.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

#define EPSILON HAMERS_EPSILON

FlowModelSubgridScaleModelSingleSpecies::FlowModelSubgridScaleModelSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_species,
    const HAMERS_SHARED_PTR<tbox::Database>& subgrid_scale_model_db):
        FlowModelSubgridScaleModel(
            object_name,
            dim,
            grid_geometry,
            num_species,
            2 + dim.getValue(),
            subgrid_scale_model_db)
{
    if (d_subgrid_scale_model_type == SUBGRID_SCALE_MODEL::VREMAN)
    {
        d_constant_sgs = Real(0.08);
        d_species_Pr_t = Real(0.9);
        
        d_constant_sgs = subgrid_scale_model_db->getRealWithDefault("constant_sgs",   d_constant_sgs);
        d_constant_sgs = subgrid_scale_model_db->getRealWithDefault("d_constant_sgs", d_constant_sgs);
        
        d_species_Pr_t = subgrid_scale_model_db->getRealWithDefault("species_Pr_t",   d_species_Pr_t);
        d_species_Pr_t = subgrid_scale_model_db->getRealWithDefault("d_species_Pr_t", d_species_Pr_t);
        
        if (subgrid_scale_model_db->keyExists("species_c_p"))
        {
            d_species_c_p = subgrid_scale_model_db->getFloat("species_c_p");
        }
        else if (subgrid_scale_model_db->keyExists("d_species_c_p"))
        {
            d_species_c_p = subgrid_scale_model_db->getFloat("d_species_c_p");
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSubgridScaleModelSingleSpecies::FlowModelSubgridScaleModelSingleSpecies()\n"
                << "No key 'species_c_p'/'d_species_c_p' found in data for subgrid scale model."
                << std::endl);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::FlowModelSubgridScaleModelSingleSpecies()\n"
            << "Unknown/unsupported subgrid scale model."
            << std::endl);
    }
}


/*
 * Return names of different derived variables required to register.
 */
std::vector<std::string>
FlowModelSubgridScaleModelSingleSpecies::getDerivedVariablesToRegister() const
{
    std::vector<std::string> var_to_register;
    var_to_register.reserve(2);
    
    var_to_register.push_back("DENSITY");
    var_to_register.push_back("VELOCITY");
    
    return var_to_register;
}


/*
 * Return different derived variables required for interpolation to sides.
 */
void
FlowModelSubgridScaleModelSingleSpecies::getDerivedVariablesForInterpolationToSideData(
    std::vector<std::string>& var_to_interpolate,
    std::vector<int>& var_to_interpolate_component_idx) const
{
    var_to_interpolate.resize(1);
    var_to_interpolate_component_idx.resize(1);
    
    var_to_interpolate[0] = "DENSITY";
    var_to_interpolate_component_idx[0] = 0;
}


/*
 * Get the variables for the derivatives used at computing subgrid scale diffusivity/viscosity at sides.
 */
void
FlowModelSubgridScaleModelSingleSpecies::getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity(
    std::vector<std::string>& derivative_var_data_str,
    std::vector<int>& derivative_var_component_idx,
    const DIRECTION::TYPE& side_direction,
    const DIRECTION::TYPE& derivative_direction)
{
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity() "
            << "not implemented for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "getCellDataOfVariablesForSideDerivativeForSubgridScaleViscosity() "
            << "not implemented for two-dimensional problem."
            << std::endl);
    }
    
    NULL_USE(side_direction);
    NULL_USE(derivative_direction);
    
    derivative_var_data_str.resize(3);
    derivative_var_component_idx.resize(3);
    
    derivative_var_data_str[0] = "VELOCITY";
    derivative_var_component_idx[0] = 0;
    
    derivative_var_data_str[1] = "VELOCITY";
    derivative_var_component_idx[1] = 1;
    
    derivative_var_data_str[2] = "VELOCITY";
    derivative_var_component_idx[2] = 2;
}


/*
 * Modify the side data of the diffusivities/viscosities at sides with subgrid scale diffusivity/viscosity.
 */
void
FlowModelSubgridScaleModelSingleSpecies::updateSideDataOfDiffusiveFluxDiffusivities(
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& var_data_for_diffusivities,
    const std::map<DIRECTION::TYPE, std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& derivatives,
    const DIRECTION::TYPE& side_direction,
    const hier::Patch& patch)
{
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "updateSideDataOfDiffusiveFluxDiffusivities()\n"
            << "not implemented for one-dimensional problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "updateSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Not implemented for two-dimensional problem."
            << std::endl);
    }
    
    if (derivatives.find(DIRECTION::X_DIRECTION) == derivatives.end())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "updateSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Cannot find derivatives in x-direction in the map of derivatives."
            << std::endl);
    }
    
    if (derivatives.find(DIRECTION::Y_DIRECTION) == derivatives.end())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "updateSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Cannot find derivatives in y-direction in the map of derivatives."
            << std::endl);
    }
    
    if (derivatives.find(DIRECTION::Z_DIRECTION) == derivatives.end())
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::"
            << "updateSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Cannot find derivatives in z-direction in the map of derivatives."
            << std::endl);
    }
    
    /*
     * Get the derivatives.
     */
    
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_x = derivatives.find(DIRECTION::X_DIRECTION)->second;
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_y = derivatives.find(DIRECTION::Y_DIRECTION)->second;
    std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > derivatives_z = derivatives.find(DIRECTION::Z_DIRECTION)->second;
    
    /*
     * Get the dimension of the interior box.
     */
    
    const hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    const hier::IntVector num_ghosts_diffus = var_data_for_diffusivities[0]->getGhostCellWidth();
    hier::IntVector ghostcell_dims_diffus = var_data_for_diffusivities[0]->getGhostBox().numberCells();
    
    const hier::IntVector num_ghosts_der = derivatives_x[0]->getGhostCellWidth();
    hier::IntVector ghostcell_dims_der = derivatives_x[0]->getGhostBox().numberCells();
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    std::vector<std::string> var_to_interpolate;
    std::vector<int> var_to_interpolate_component_idx;
    getDerivedVariablesForInterpolationToSideData(var_to_interpolate, var_to_interpolate_component_idx);
    
    TBOX_ASSERT(static_cast<int>(var_data_for_diffusivities.size()) == 3 + d_dim.getValue() +
        static_cast<int>(var_to_interpolate));
    
    for (int vi = 0; vi < static_cast<int>(var_data_for_diffusivities.size()); vi++)
    {
        TBOX_ASSERT(var_data_for_diffusivities[vi]->getGhostCellWidth() == num_ghosts_diffus);
        TBOX_ASSERT(var_data_for_diffusivities[vi]->getGhostBox().contains(interior_box));
    }
    
    TBOX_ASSERT(static_cast<int>(derivatives_x.size()) == 3);
    TBOX_ASSERT(static_cast<int>(derivatives_y.size()) == 3);
    TBOX_ASSERT(static_cast<int>(derivatives_z.size()) == 3);
    
    for (int vi = 0; vi < static_cast<int>(derivatives_x.size()); vi++)
    {
        TBOX_ASSERT(derivatives_x[vi]->getGhostCellWidth() == num_ghosts_der);
        TBOX_ASSERT(derivatives_x[vi]->getGhostBox().contains(interior_box));
        
        TBOX_ASSERT(derivatives_y[vi]->getGhostCellWidth() == num_ghosts_der);
        TBOX_ASSERT(derivatives_y[vi]->getGhostBox().contains(interior_box));
        
        TBOX_ASSERT(derivatives_z[vi]->getGhostCellWidth() == num_ghosts_der);
        TBOX_ASSERT(derivatives_z[vi]->getGhostBox().contains(interior_box));
    }
#endif
    
    // Get the grid spacing.
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    // Only Dimension(3) is implemented.
    
    int side_normal = 0;
    switch(side_direction)
    {
        case DIRECTION::X_DIRECTION:
        {
            side_normal = 0;
            break;
        }
        case DIRECTION::Y_DIRECTION:
        {
            
            side_normal = 1;
            break;
        }
        case DIRECTION::Z_DIRECTION:
        {
            side_normal = 2;
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelSubgridScaleModelSingleSpecies::"
                << "updateSideDataOfDiffusiveFluxDiffusivities()\n"
                << "There are only x-direction, y-direction and z-direction for three-dimensional problem."
                << std::endl);
        }
    }
    
    hier::IntVector domain_lo(d_dim);
    hier::IntVector domain_dims(d_dim);
    
    ghostcell_dims_diffus[side_normal]++;
    ghostcell_dims_der[side_normal]++;
    
    const Real delta = Real(dx[side_normal]);
    
    domain_lo = hier::IntVector::getZero(d_dim);
    domain_dims = interior_dims;
    
    domain_lo[side_normal] -= num_ghosts_diffus[side_normal];
    domain_dims[side_normal] += 1;
    domain_dims[side_normal] += 2*(num_ghosts_diffus[side_normal]);
    
    Real* mu    = var_data_for_diffusivities[0]->getPointer(side_normal, 0);
    // Real* mu_v  = var_data_for_diffusivities[1]->getPointer(side_normal, 0);
    Real* kappa = var_data_for_diffusivities[2]->getPointer(side_normal, 0);
    
    const Real* const rho = var_data_for_diffusivities[6 + 0]->getPointer(side_normal, 0);
    
    const Real* const ddx_u = derivatives_x[0]->getPointer(side_normal, 0);
    const Real* const ddx_v = derivatives_x[1]->getPointer(side_normal, 0);
    const Real* const ddx_w = derivatives_x[2]->getPointer(side_normal, 0);
    const Real* const ddy_u = derivatives_y[0]->getPointer(side_normal, 0);
    const Real* const ddy_v = derivatives_y[1]->getPointer(side_normal, 0);
    const Real* const ddy_w = derivatives_y[2]->getPointer(side_normal, 0);
    const Real* const ddz_u = derivatives_z[0]->getPointer(side_normal, 0);
    const Real* const ddz_v = derivatives_z[1]->getPointer(side_normal, 0);
    const Real* const ddz_w = derivatives_z[2]->getPointer(side_normal, 0);
    
    if (d_subgrid_scale_model_type == SUBGRID_SCALE_MODEL::VREMAN)
    {
        updateSideDataOfDiffusiveFluxDiffusivitiesVreman(
            mu,
            kappa,
            rho,
            ddx_u,
            ddx_v,
            ddx_w,
            ddy_u,
            ddy_v,
            ddy_w,
            ddz_u,
            ddz_v,
            ddz_w,
            num_ghosts_diffus,
            num_ghosts_der,
            ghostcell_dims_diffus,
            ghostcell_dims_der,
            domain_lo,
            domain_dims,
            delta);
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::updateSideDataOfDiffusiveFluxDiffusivities()\n"
            << "Unknown/unsupported subgrid scale model."
            << std::endl);
    }
}


/*
 * Put the characteristics of this class into the restart database.
 */
void
FlowModelSubgridScaleModelSingleSpecies::putToRestart(
    const HAMERS_SHARED_PTR<tbox::Database>& restart_subgrid_scale_model_db) const
{
    putToRestartBase(restart_subgrid_scale_model_db);
    
    if (d_subgrid_scale_model_type == SUBGRID_SCALE_MODEL::VREMAN)
    {
        restart_subgrid_scale_model_db->putDouble("d_constant_sgs", d_constant_sgs);
        restart_subgrid_scale_model_db->putDouble("d_species_Pr_t", d_species_Pr_t);
        restart_subgrid_scale_model_db->putDouble("d_species_c_p",  d_species_c_p);
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelSubgridScaleModelSingleSpecies::putToRestart()\n"
            << "Unknown/unsupported subgrid scale model."
            << std::endl);
    }
}


/*
 * Kernal to modify the side data of the diffusivities/viscosities at sides with Vreman's subgrid scale
 * diffusivity/viscosity.
 */
void
FlowModelSubgridScaleModelSingleSpecies::updateSideDataOfDiffusiveFluxDiffusivitiesVreman(
    Real* mu,
    Real* kappa,
    const Real* const rho,
    const Real* const ddx_u,
    const Real* const ddx_v,
    const Real* const ddx_w,
    const Real* const ddy_u,
    const Real* const ddy_v,
    const Real* const ddy_w,
    const Real* const ddz_u,
    const Real* const ddz_v,
    const Real* const ddz_w,
    const hier::IntVector& num_ghosts_diffus,
    const hier::IntVector& num_ghosts_der,
    const hier::IntVector& ghostcell_dims_diffus,
    const hier::IntVector& ghostcell_dims_der,
    const hier::IntVector& domain_lo,
    const hier::IntVector& domain_dims,
    const Real& delta) const
{
    /*
     * Get the local lower index, numbers of cells in each dimension and number of ghost cells.
     */
    
    const int domain_lo_0 = domain_lo[0];
    const int domain_lo_1 = domain_lo[1];
    const int domain_lo_2 = domain_lo[2];
    const int domain_dim_0 = domain_dims[0];
    const int domain_dim_1 = domain_dims[1];
    const int domain_dim_2 = domain_dims[2];
    
    const int num_ghosts_0_diffus = num_ghosts_diffus[0];
    const int num_ghosts_1_diffus = num_ghosts_diffus[1];
    const int num_ghosts_2_diffus = num_ghosts_diffus[2];
    const int ghostcell_dim_0_diffus = ghostcell_dims_diffus[0];
    const int ghostcell_dim_1_diffus = ghostcell_dims_diffus[1];
    
    const int num_ghosts_0_der = num_ghosts_der[0];
    const int num_ghosts_1_der = num_ghosts_der[1];
    const int num_ghosts_2_der = num_ghosts_der[2];
    const int ghostcell_dim_0_der = ghostcell_dims_der[0];
    const int ghostcell_dim_1_der = ghostcell_dims_der[1];
    
    for (int k = domain_lo_2; k < domain_lo_2 + domain_dim_2; k++)
    {
        for (int j = domain_lo_1; j < domain_lo_1 + domain_dim_1; j++)
        {
            HAMERS_PRAGMA_SIMD
            for (int i = domain_lo_0; i < domain_lo_0 + domain_dim_0; i++)
            {
                // Compute the linear indices.
                const int idx_diffus = (i + num_ghosts_0_diffus) +
                    (j + num_ghosts_1_diffus)*ghostcell_dim_0_diffus +
                    (k + num_ghosts_2_diffus)*ghostcell_dim_0_diffus*
                        ghostcell_dim_1_diffus;
                
                const int idx_der = (i + num_ghosts_0_der) +
                    (j + num_ghosts_1_der)*ghostcell_dim_0_der +
                    (k + num_ghosts_2_der)*ghostcell_dim_0_der*
                        ghostcell_dim_1_der;
                
                const Real g_11 = ddx_u[idx_der]*ddx_u[idx_der] + ddy_u[idx_der]*ddy_u[idx_der] + ddz_u[idx_der]*ddz_u[idx_der];
                const Real g_22 = ddx_v[idx_der]*ddx_v[idx_der] + ddy_v[idx_der]*ddy_v[idx_der] + ddz_v[idx_der]*ddz_v[idx_der];
                const Real g_33 = ddx_w[idx_der]*ddx_w[idx_der] + ddy_w[idx_der]*ddy_w[idx_der] + ddz_w[idx_der]*ddz_w[idx_der];
                
                const Real g_12 = ddx_u[idx_der]*ddx_v[idx_der] + ddy_u[idx_der]*ddy_v[idx_der] + ddz_u[idx_der]*ddz_v[idx_der];
                const Real g_13 = ddx_u[idx_der]*ddx_w[idx_der] + ddy_u[idx_der]*ddy_w[idx_der] + ddz_u[idx_der]*ddz_w[idx_der];
                const Real g_23 = ddx_v[idx_der]*ddx_w[idx_der] + ddy_v[idx_der]*ddy_w[idx_der] + ddz_v[idx_der]*ddz_w[idx_der];
                
                Real sig_D = g_11*g_22 - g_12*g_12 + g_11*g_33 - g_13*g_13 + g_22*g_33 - g_23*g_23;
                sig_D = std::max(sig_D, Real(0));
                sig_D /= (
                    ddx_u[idx_der]*ddx_u[idx_der] + ddy_u[idx_der]*ddy_u[idx_der] + ddz_u[idx_der]*ddz_u[idx_der] +
                    ddx_v[idx_der]*ddx_v[idx_der] + ddy_v[idx_der]*ddy_v[idx_der] + ddz_v[idx_der]*ddz_v[idx_der] +
                    ddx_w[idx_der]*ddx_w[idx_der] + ddy_w[idx_der]*ddy_w[idx_der] + ddz_w[idx_der]*ddz_w[idx_der] +
                    EPSILON);
                
                sig_D = d_constant_sgs*std::sqrt(sig_D + EPSILON);
                
                const Real mu_SGS    = sig_D*delta*delta*rho[idx_diffus];
                const Real kappa_SGS = d_species_c_p*mu_SGS/d_species_Pr_t;
                
                mu[idx_diffus]    += mu_SGS;
                kappa[idx_diffus] += kappa_SGS;
            }
        }
    }
}
