#ifndef REFINE_REGIONS_TAGGER_HPP
#define REFINE_REGIONS_TAGGER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "util/refinement_regions_taggers/RefineBox.hpp"
#include "util/refinement_regions_taggers/RefineCircle.hpp"
#include "util/refinement_regions_taggers/RefineTriangle.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/pdat/CellData.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class RefineRegionsTagger
{
    public:
        RefineRegionsTagger(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const HAMERS_SHARED_PTR<tbox::Database>& refine_regions_tagger_db,
            const bool is_db_from_restart,
            const bool read_on_restart);
        
        /*
         * Put the characteristics of the refine regions tagger class into the restart database.
         */
        void
        putToRestart(const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Get the characteristics of the refine regions tagger class from the restart database.
         */
        void
        getFromRestart(const HAMERS_SHARED_PTR<tbox::Database>& restart_db);
        
        /*
         * Tag cells on a patch for refinement.
         */
        void
        tagCellsOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags) const;
        
    private:
        /*
         * Tag cells on a patch for refinement using refine boxes.
         */
        void
        tagCellsOnPatchUsingRefineBoxes(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags) const;
        
        /*
         * Tag cells on a patch for refinement using refine circles.
         */
        void
        tagCellsOnPatchUsingRefineCircles(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags) const;
        
        /*
         * Tag cells on a patch for refinement using refine triangles.
         */
        void
        tagCellsOnPatchUsingRefineTriangles(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags) const;
        
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
         * Refine boxes settings.
         */
        bool d_use_refine_boxes;
        int d_num_refine_boxes;
        std::vector<RefineBox> d_refine_boxes;
        
        /*
         * Refine circles settings.
         */
        bool d_use_refine_circles;
        int d_num_refine_circles;
        std::vector<RefineCircle> d_refine_circles;
        
        /*
         * Refine triangles settings.
         */
        bool d_use_refine_triangles;
        int d_num_refine_triangles;
        std::vector<RefineTriangle> d_refine_triangles;
};

#endif /* REFINE_REGIONS_TAGGER_HPP */
