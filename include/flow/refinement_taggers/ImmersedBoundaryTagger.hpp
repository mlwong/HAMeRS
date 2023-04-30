#ifndef IMMERSED_BOUNDARY_TAGGER_HPP
#define IMMERSED_BOUNDARY_TAGGER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Dimension.h"

using namespace SAMRAI;

class ImmersedBoundaryTagger
{
    public:
        ImmersedBoundaryTagger(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const hier::IntVector& num_cells_buffer);
        
        /*
         * Get the number of ghost cells needed by the immersed boundary tagger.
         */
        hier::IntVector
        getImmersedBoundaryTaggerNumberOfGhostCells() const
        {
            return d_num_immersed_boundary_tagger_ghosts;
        }
        
        /*
         * Print all characteristics of the immersed boundary tagger class.
         */
        void
        printClassData(std::ostream& os) const;
        
        // /*
        //  * Put the characteristics of the immersed boundary tagger class into the restart database.
        //  */
        // void
        // putToRestart(
        //     const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Tag cells on a patch for refinement using gradient sensors.
         */
        void
        tagCellsOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
    private:
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
         * Number of cells in buffer zone from immersed boundaries for tagging.
         */
        hier::IntVector d_num_cells_buffer;
        
        /*
         * Number of ghost cells needed by the immersed boundary tagger.
         */
        hier::IntVector d_num_immersed_boundary_tagger_ghosts;
        
};


#endif /* IMMERSED_BOUNDARY_TAGGER_HPP */
