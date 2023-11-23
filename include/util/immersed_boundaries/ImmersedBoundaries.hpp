#ifndef IMMERSED_BOUNDARIES_HPP
#define IMMERSED_BOUNDARIES_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/pdat/CellData.h"

#include <string>

#define IB_EPSILON HAMERS_EPSILON

using namespace SAMRAI;

namespace IB_MASK
{
    enum TYPE { FLUID = 0,
                IB_GHOST,
                BODY };
}

class ImmersedBoundaries
{
    public:
        ImmersedBoundaries(
            const std::string& object_name,
            const std::string& project_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<tbox::Database>& initial_conditions_db,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry):
                d_object_name(object_name),
                d_project_name(project_name),
                d_dim(dim),
                d_initial_conditions_db(initial_conditions_db),
                d_grid_geometry(grid_geometry),
                d_num_immersed_boundary_ghosts(-hier::IntVector::getOne(dim))
                
        {}
        
        virtual ~ImmersedBoundaries() {}
        
        void setNumberOfImmersedBoundaryGhosts(const hier::IntVector& num_immersed_boundary_ghosts)
        {
            d_num_immersed_boundary_ghosts = num_immersed_boundary_ghosts;
            
            if (d_dim == tbox::Dimension(2))
            {
                if (d_num_immersed_boundary_ghosts[0] != d_num_immersed_boundary_ghosts[1])
                {
                    TBOX_ERROR(d_object_name
                        << ": ImmersedBoundaries::"
                        << "setNumberOfImmersedBoundaryGhosts()\n"
                        << "Can only have same numbers of immersed boundary ghosts in all directions."
                        << std::endl);
                }
            }
            else if (d_dim == tbox::Dimension(3))
            {
                if ((d_num_immersed_boundary_ghosts[0] != d_num_immersed_boundary_ghosts[1]) ||
                    (d_num_immersed_boundary_ghosts[1] != d_num_immersed_boundary_ghosts[2]))
                {
                    TBOX_ERROR(d_object_name
                        << ": ImmersedBoundaries::"
                        << "setNumberOfImmersedBoundaryGhosts()\n"
                        << "Can only have same numbers of immersed boundary ghosts in all directions."
                        << std::endl);
                }
            }
        }
        
        void setImmersedBoundaryVariablesOnPatch(
            const hier::Patch& patch,
            const double data_time,
            const bool initial_time,
            const hier::Box& domain,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_mask,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_wall_distance,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_surface_normal)
        {
            NULL_USE(data_time);
            
            const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
                HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                    patch.getPatchGeometry()));
            
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(patch_geom);
#endif
            
            // Get the dimensions of box that covers the interior of patch.
            hier::Box interior_box = patch.getBox();
            
            const hier::IntVector num_ghosts = data_mask->getGhostCellWidth();
            const hier::IntVector ghostcell_dims = data_mask->getGhostBox().numberCells();
            
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(num_ghosts == data_wall_distance->getGhostCellWidth());
            TBOX_ASSERT(num_ghosts == data_surface_normal->getGhostCellWidth());
#endif
            
            /*
             * Get the local lower index and number of cells in each direction of the domain.
             */
            
            hier::IntVector domain_lo(d_dim);
            hier::IntVector domain_dims(d_dim);
            
            if (domain.empty())
            {
                domain_lo = -num_ghosts;
                domain_dims = ghostcell_dims;
            }
            else
            {
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(data_mask.contains(domain));
                TBOX_ASSERT(data_wall_distance.contains(domain));
                TBOX_ASSERT(data_surface_normal.contains(domain));
#endif
                
                domain_lo = domain.lower() - interior_box.lower();
                domain_dims = domain.numberCells();
            }
            
            setImmersedBoundaryVariablesOnPatch(
                patch,
                data_time,
                initial_time,
                domain_lo,
                domain_dims,
                data_mask,
                data_wall_distance,
                data_surface_normal);
        }
        
        void setImmersedBoundaryVariablesOnPatch(
            const hier::Patch& patch,
            const double data_time,
            const bool initial_time,
            const hier::IntVector& domain_lo,
            const hier::IntVector& domain_dims,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& data_mask,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_wall_distance,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_surface_normal);
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Name of the project.
         */
        std::string d_project_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Initial conditions database.
         */
        const HAMERS_SHARED_PTR<tbox::Database> d_initial_conditions_db;

        /*
         * Number of immersed boundary ghost cells (IB_MASK::IBGHOST).
         */
        hier::IntVector d_num_immersed_boundary_ghosts;
        
};

#endif /* IMMERSED_BOUNDARIES_HPP */
