#include "util/refinement_regions_taggers/RefineRegionsTagger.hpp"

// Compute the area of a triangle.
static inline __attribute__((always_inline)) Real computeTriangleArea(
    const Real x_0, const Real y_0,
    const Real x_1, const Real y_1,
    const Real x_2, const Real y_2)
{
    return std::abs(Real(1)/Real(2)*(x_0*(y_1 - y_2) + x_1*(y_2 - y_0) + x_2*(y_0 - y_1)));
}


RefineRegionsTagger::RefineRegionsTagger(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const HAMERS_SHARED_PTR<tbox::Database>& refine_regions_tagger_db,
    const bool is_db_from_restart):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry)
{
    if (refine_regions_tagger_db != nullptr)
    {
        if (!is_db_from_restart)
        {
            /*
             * Get the refine boxes in the sequence of 'refine_box_0', 'refine_box_1', ...
             */
            
            d_use_refine_boxes = false;
            int bi = 0;
            bool refine_box_exists = false;
            std::string refine_box_name = "refine_box_" + std::to_string(bi);
            refine_box_exists = refine_regions_tagger_db->keyExists(refine_box_name);
            
            if (refine_box_exists)
            {
                d_use_refine_boxes = true;
            }
            
            while (refine_box_exists)
            {
                std::shared_ptr<tbox::Database> refine_box_db = refine_regions_tagger_db->getDatabase(refine_box_name);
                
                // Get the lower and upper coordinates of the refine box.
                std::vector<Real> x_lo;
                std::vector<Real> x_hi;
                
                if (refine_box_db->keyExists("x_lo"))
                {
                    x_lo = refine_box_db->getRealVector("x_lo");
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Key data 'x_lo' not found for '" << refine_box_name << "' in input database."
                        << std::endl);
                }
                
                if (refine_box_db->keyExists("x_hi"))
                {
                    x_hi = refine_box_db->getRealVector("x_hi");
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Key data 'x_hi' not found for '" << refine_box_name << "' in input database."
                        << std::endl);
                }
                
                // Check the size of the vectors.
                
                if (static_cast<int>(x_lo.size()) != d_dim.getValue())
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "The size of the vector 'x_lo' for '" << refine_box_name << "' is not correct."
                        << std::endl);
                }
                
                if (static_cast<int>(x_hi.size()) != d_dim.getValue())
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "The size of the vector 'x_hi' for '" << refine_box_name << "' is not correct."
                        << std::endl);
                }
                
                // Check that the lower coordinates are less than or equal to the upper coordinates.
                for (int i = 0; i < d_dim.getValue(); i++)
                {
                    if (x_lo[i] > x_hi[i])
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "The lower coordinates are greater than the upper coordinates for '" << refine_box_name << "'."
                            << std::endl);
                    }
                }
                
                // Get the number of refine levels.
                int num_refine_levels = 0;
                
                if (refine_box_db->keyExists("num_refine_levels"))
                {
                    num_refine_levels = refine_box_db->getInteger("num_refine_levels");
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Key data 'num_refine_levels' not found for '" << refine_box_name << "' in input database."
                        << std::endl);
                }
                
                // Check that the number of refine levels is positive.
                if (num_refine_levels <= 0)
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "The number of refine levels is not positive for '" << refine_box_name << "'."
                        << std::endl);
                }
                
                // Add the box to the refine boxes.
                d_refine_boxes.push_back(RefineBox(x_lo, x_hi, num_refine_levels));
                
                // Search for next box.
                bi++;
                refine_box_name = "refine_box_" + std::to_string(bi);
                refine_box_exists = refine_regions_tagger_db->keyExists(refine_box_name);
            }
            
            d_num_refine_boxes = static_cast<int>(d_refine_boxes.size());
            
            /*
             * Get the refine triangles in the sequence of 'refine_triangle_0', 'refine_triangle_1', ...
             */
            
            d_use_refine_triangles = false;
            int ti = 0;
            bool refine_triangle_exists = false;
            std::string refine_triangle_name = "refine_triangle_" + std::to_string(ti);
            refine_triangle_exists = refine_regions_tagger_db->keyExists(refine_triangle_name);
            
            if (refine_triangle_exists)
            {
                d_use_refine_triangles = true;
            }
            
            // Check whether the problem is one-dimensional.
            if (d_dim == tbox::Dimension(1))
            {
                TBOX_ERROR(d_object_name
                    << ": "
                    << "Refine triangles are not allowed in one-dimensional problem."
                    << std::endl);
            }
            
            while (refine_triangle_exists)
            {
                std::shared_ptr<tbox::Database> refine_triangle_db = refine_regions_tagger_db->getDatabase(refine_triangle_name);
                
                // Get the coordinates of the refine triangle.
                std::vector<Real> point_coord_0;
                std::vector<Real> point_coord_1;
                std::vector<Real> point_coord_2;
                
                if (refine_triangle_db->keyExists("point_coord_0"))
                {
                    point_coord_0 = refine_triangle_db->getRealVector("point_coord_0");
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Key data 'point_coord_0' not found for '" << refine_triangle_name << "' in input database."
                        << std::endl);
                }
                
                if (refine_triangle_db->keyExists("point_coord_1"))
                {
                    point_coord_1 = refine_triangle_db->getRealVector("point_coord_1");
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Key data 'point_coord_1' not found for '" << refine_triangle_name << "' in input database."
                        << std::endl);
                }
                
                if (refine_triangle_db->keyExists("point_coord_2"))
                {
                    point_coord_2 = refine_triangle_db->getRealVector("point_coord_2");
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Key data 'point_coord_2' not found for '" << refine_triangle_name << "' in input database."
                        << std::endl);
                }
                
                // Check that the size of the vectors is correct.
                if (static_cast<int>(point_coord_0.size()) != d_dim.getValue())
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "The size of the vector 'point_coord_0' for '" << refine_triangle_name << "' is not correct."
                        << std::endl);
                }
                
                if (static_cast<int>(point_coord_1.size()) != d_dim.getValue())
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "The size of the vector 'point_coord_1' for '" << refine_triangle_name << "' is not correct."
                        << std::endl);
                }
                
                if (static_cast<int>(point_coord_2.size()) != d_dim.getValue())
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "The size of the vector 'point_coord_2' for '" << refine_triangle_name << "' is not correct."
                        << std::endl);
                }
                
                // Get the number of refine levels.
                int num_refine_levels = 0;
                
                if (refine_triangle_db->keyExists("num_refine_levels"))
                {
                    num_refine_levels = refine_triangle_db->getInteger("num_refine_levels");
                }
                else
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "Key data 'num_refine_levels' not found for '" << refine_triangle_name << "' in input database."
                        << std::endl);
                }
                
                // Get the direction if the problem is three-dimensional.
                
                DIRECTION::TYPE direction = DIRECTION::X_DIRECTION;
                if (d_dim == tbox::Dimension(3))
                {
                    if (refine_triangle_db->keyExists("direction"))
                    {
                        const std::string direction_str = refine_triangle_db->getString("direction");
                        if (direction_str == "x" || direction_str == "X")
                        {
                            direction = DIRECTION::X_DIRECTION;
                        }
                        else if (direction_str == "y" || direction_str == "Y")
                        {
                            direction = DIRECTION::Y_DIRECTION;
                        }
                        else if (direction_str == "z" || direction_str == "Z")
                        {
                            direction = DIRECTION::Z_DIRECTION;
                        }
                        else
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "Unknown direction '" << direction_str << "' for '" << refine_triangle_name << "'."
                                << std::endl);
                        }
                    }
                    else
                    {
                        TBOX_ERROR(d_object_name
                            << ": "
                            << "Key data 'direction' not found for '" << refine_triangle_name << "' in input database."
                            << std::endl);
                    }
                }
                
                // Make sure the triangle is on a plane if the problem is three-dimensional.
                if (d_dim == tbox::Dimension(3))
                {
                    if (direction == DIRECTION::X_DIRECTION)
                    {
                        if (std::abs(point_coord_0[0] - point_coord_1[0]) > std::numeric_limits<Real>::epsilon() ||
                            std::abs(point_coord_0[0] - point_coord_2[0]) > std::numeric_limits<Real>::epsilon() ||
                            std::abs(point_coord_1[0] - point_coord_2[0]) > std::numeric_limits<Real>::epsilon())
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The triangle '" << refine_triangle_name << "' is not on a plane."
                                << std::endl);
                        }
                    }
                    else if (direction == DIRECTION::Y_DIRECTION)
                    {
                        if (std::abs(point_coord_0[1] - point_coord_1[1]) > std::numeric_limits<Real>::epsilon() ||
                            std::abs(point_coord_0[1] - point_coord_2[1]) > std::numeric_limits<Real>::epsilon() ||
                            std::abs(point_coord_1[1] - point_coord_2[1]) > std::numeric_limits<Real>::epsilon())
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The triangle '" << refine_triangle_name << "' is not on a plane."
                                << std::endl);
                        }
                    }
                    else if (direction == DIRECTION::Z_DIRECTION)
                    {
                        if (std::abs(point_coord_0[2] - point_coord_1[2]) > std::numeric_limits<Real>::epsilon() ||
                            std::abs(point_coord_0[2] - point_coord_2[2]) > std::numeric_limits<Real>::epsilon() ||
                            std::abs(point_coord_1[2] - point_coord_2[2]) > std::numeric_limits<Real>::epsilon())
                        {
                            TBOX_ERROR(d_object_name
                                << ": "
                                << "The triangle '" << refine_triangle_name << "' is not on a plane."
                                << std::endl);
                        }
                    }
                }
                
                // Check the number of refine levels is positive.
                if (num_refine_levels <= 0)
                {
                    TBOX_ERROR(d_object_name
                        << ": "
                        << "The number of refine levels is not positive for '" << refine_triangle_name << "'."
                        << std::endl);
                }
                
                // Add the triangle to the refine triangles.
                d_refine_triangles.push_back(RefineTriangle(point_coord_0, point_coord_1, point_coord_2, direction, num_refine_levels));
                
                // Search for next refine trapezoid.
                ti++;
                refine_triangle_name = "refine_triangle_" + std::to_string(ti);
                refine_triangle_exists = refine_regions_tagger_db->keyExists(refine_triangle_name);
            }
            
            d_num_refine_triangles = static_cast<int>(d_refine_triangles.size());
        }
        else
        {
            getFromRestart(refine_regions_tagger_db);
        }
    }
    else
    {
        TBOX_WARNING(d_object_name
            << ": "
            << "Key data 'Refine_regions_tagger' not found in input/restart database."
            << " No refinement with user-defined refine regions will occur."
            << std::endl);
    }
}


/*
 * Put the characteristics of the refine regions tagger class into the restart database.
 */
void
RefineRegionsTagger::putToRestart(const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
{
    TBOX_ASSERT(restart_db);
    
    /*
     * Put the refine boxes into the restart database.
     */
    
    restart_db->putBool("d_use_refine_boxes", d_use_refine_boxes);
    restart_db->putInteger("d_num_refine_boxes", d_num_refine_boxes);
    for (int ib = 0; ib < static_cast<int>(d_refine_boxes.size()); ib++)
    {
        std::string refine_box_name = "refine_box_" + std::to_string(ib);
        std::shared_ptr<tbox::Database> refine_box_db = restart_db->putDatabase(refine_box_name);
        
        refine_box_db->putRealVector("x_lo", d_refine_boxes[ib].x_lo);
        refine_box_db->putRealVector("x_hi", d_refine_boxes[ib].x_hi);
        refine_box_db->putInteger("num_refine_levels", d_refine_boxes[ib].num_refine_levels);
    }
    
    /*
     * Put the refine triangles into the restart database.
     */
    
    restart_db->putBool("d_use_refine_triangles", d_use_refine_triangles);
    restart_db->putInteger("d_num_refine_triangles", d_num_refine_triangles);
    for (int it = 0; it < static_cast<int>(d_refine_triangles.size()); it++)
    {
        std::string refine_triangle_name = "refine_triangle_" + std::to_string(it);
        std::shared_ptr<tbox::Database> refine_triangle_db = restart_db->putDatabase(refine_triangle_name);
        
        refine_triangle_db->putRealVector("point_coord_0", d_refine_triangles[it].point_coord_0);
        refine_triangle_db->putRealVector("point_coord_1", d_refine_triangles[it].point_coord_1);
        refine_triangle_db->putRealVector("point_coord_2", d_refine_triangles[it].point_coord_2);
        
        const DIRECTION::TYPE direction = d_refine_triangles[it].direction;
        if (direction == DIRECTION::X_DIRECTION)
        {
            refine_triangle_db->putString("direction", "x");
        }
        else if (direction == DIRECTION::Y_DIRECTION)
        {
            refine_triangle_db->putString("direction", "y");
        }
        else if (direction == DIRECTION::Z_DIRECTION)
        {
            refine_triangle_db->putString("direction", "z");
        }
        refine_triangle_db->putInteger("num_refine_levels", d_refine_triangles[it].num_refine_levels);
    }
}


/*
 * Get the characteristics of the refine regions tagger class from the restart database.
 */
void
RefineRegionsTagger::getFromRestart(const HAMERS_SHARED_PTR<tbox::Database>& restart_db)
{
    TBOX_ASSERT(restart_db);
    
    /*
     * Get the refine boxes from the restart database.
     */
    
    d_use_refine_boxes = restart_db->getBool("d_use_refine_boxes");
    d_num_refine_boxes = restart_db->getInteger("d_num_refine_boxes");
    for (int ib = 0; ib < d_num_refine_boxes; ib++)
    {
        std::string refine_box_name = "refine_box_" + std::to_string(ib);
        std::shared_ptr<tbox::Database> refine_box_db = restart_db->getDatabase(refine_box_name);
        
        // Get the lower and upper coordinates of the refine box.
        std::vector<Real> x_lo;
        std::vector<Real> x_hi;
        x_lo = refine_box_db->getRealVector("x_lo");
        x_hi = refine_box_db->getRealVector("x_hi");
        const int num_refine_levels = refine_box_db->getInteger("num_refine_levels");
        
        d_refine_boxes.push_back(RefineBox(x_lo, x_hi, num_refine_levels));
    }
    
    /*
     * Get the refine triangles from the restart database.
     */
    
    d_use_refine_triangles = restart_db->getBool("d_use_refine_triangles");
    d_num_refine_triangles = restart_db->getInteger("d_num_refine_triangles");
    for (int it = 0; it < d_num_refine_triangles; it++)
    {
        std::string refine_triangle_name = "refine_triangle_" + std::to_string(it);
        std::shared_ptr<tbox::Database> refine_triangle_db = restart_db->getDatabase(refine_triangle_name);
        
        // Get the coordinates of the three points of the refine triangle.
        std::vector<Real> point_coord_0;
        std::vector<Real> point_coord_1;
        std::vector<Real> point_coord_2;
        point_coord_0 = refine_triangle_db->getRealVector("point_coord_0");
        point_coord_1 = refine_triangle_db->getRealVector("point_coord_1");
        point_coord_2 = refine_triangle_db->getRealVector("point_coord_2");
        
        // Get the direction of the refine triangle.
        const std::string direction_str = refine_triangle_db->getString("direction");
        DIRECTION::TYPE direction = DIRECTION::Z_DIRECTION;
        if (direction_str == "x")
        {
            direction = DIRECTION::X_DIRECTION;
        }
        else if (direction_str == "y")
        {
            direction = DIRECTION::Y_DIRECTION;
        }
        else if (direction_str == "z")
        {
            direction = DIRECTION::Z_DIRECTION;
        }
        else
        {
            TBOX_ERROR(d_object_name
                << ": RefineRegionsTagger::getFromRestart()\n"
                << "Unknown direction of refine triangle encountered in input file."
                << std::endl);
        }
        
        const int num_refine_levels = refine_triangle_db->getInteger("num_refine_levels");
        
        d_refine_triangles.push_back(RefineTriangle(point_coord_0, point_coord_1, point_coord_2, direction, num_refine_levels));
    }
}


/*
 * Tag cells on a patch for refinement.
 */
void
RefineRegionsTagger::tagCellsOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags) const
{
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(tags->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    if (d_use_refine_boxes)
    {
        tagCellsOnPatchUsingRefineBoxes(patch, tags);
    }
    
    if (d_use_refine_triangles)
    {
        tagCellsOnPatchUsingRefineTriangles(patch, tags);
    }
}


/*
 * Tag cells on a patch for refinement using refine boxes.
 */
void
RefineRegionsTagger::tagCellsOnPatchUsingRefineBoxes(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags) const
{
    // Get the patch level and patch dimensions.
    const int patch_level =  patch.getPatchLevelNumber();
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* const dx = patch_geom->getDx();
    const double* const patch_xlo = patch_geom->getXLower();
    
    // Get the dimensions of box that covers the interior of Patch.
    hier::Box patch_box = patch.getBox();
    const hier::IntVector patch_dims = patch_box.numberCells();
    
    // Get the pointer to the tags.
    int* tag_ptr = tags->getPointer(0);
    
    const Real half = Real(1)/Real(2);
    
    if (d_dim == tbox::Dimension(1))
    {
        const int patch_dim_0 = patch_dims[0];
        
        const Real patch_xlo_0 = Real(patch_xlo[0]);
        
        const Real dx_0 = Real(dx[0]);
        
        // Loop over the refine boxes to tag cells for refinement.
        for (int ib = 0; ib < static_cast<int>(d_refine_boxes.size()); ib++)
        {
            // Get the lower and upper coordinates of the refine box.
            const std::vector<Real>& refine_box_x_lo = d_refine_boxes[ib].x_lo;
            const std::vector<Real>& refine_box_x_hi = d_refine_boxes[ib].x_hi;
            
            const Real refine_box_x_lo_0 = refine_box_x_lo[0];
            const Real refine_box_x_hi_0 = refine_box_x_hi[0];
            
            // Get the number of refine levels.
            const int num_refine_levels = d_refine_boxes[ib].num_refine_levels;
            
            if (patch_level < num_refine_levels)
            {
                HAMERS_PRAGMA_SIMD
                for (int i = 0; i < patch_dim_0; i++)
                {
                    // Compute the coordinate.
                    const Real x = patch_xlo_0 + (Real(i) + half)*dx_0;
                    
                    if (x >= refine_box_x_lo_0 && x <= refine_box_x_hi_0)
                    {
                        tag_ptr[i] = 1;
                    }
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int patch_dim_0 = patch_dims[0];
        const int patch_dim_1 = patch_dims[1];
        
        const Real patch_xlo_0 = Real(patch_xlo[0]);
        const Real patch_xlo_1 = Real(patch_xlo[1]);
        
        const Real dx_0 = Real(dx[0]);
        const Real dx_1 = Real(dx[1]);
        
        // Loop over the refine boxes to tag cells for refinement.
        for (int ib = 0; ib < d_num_refine_boxes; ib++)
        {
            // Get the lower and upper coordinates of the refine box.
            const std::vector<Real>& refine_box_x_lo = d_refine_boxes[ib].x_lo;
            const std::vector<Real>& refine_box_x_hi = d_refine_boxes[ib].x_hi;
            
            const Real refine_box_x_lo_0 = refine_box_x_lo[0];
            const Real refine_box_x_lo_1 = refine_box_x_lo[1];
            const Real refine_box_x_hi_0 = refine_box_x_hi[0];
            const Real refine_box_x_hi_1 = refine_box_x_hi[1];
            
            // Get the number of refine levels.
            const int num_refine_levels = d_refine_boxes[ib].num_refine_levels;
            
            if (patch_level < num_refine_levels)
            {
                for (int j = 0; j < patch_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < patch_dim_0; i++)
                    {
                        // Compute index into linear data array.
                        int idx_cell = i + j*patch_dim_0;
                        
                        // Compute the coordinates.
                        Real x[2];
                        x[0] = patch_xlo_0 + (Real(i) + half)*dx_0;
                        x[1] = patch_xlo_1 + (Real(j) + half)*dx_1;
                        
                        if (x[0] >= refine_box_x_lo_0 && x[0] <= refine_box_x_hi_0 &&
                            x[1] >= refine_box_x_lo_1 && x[1] <= refine_box_x_hi_1)
                        {
                            tag_ptr[idx_cell] = 1;
                        }
                    }
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int patch_dim_0 = patch_dims[0];
        const int patch_dim_1 = patch_dims[1];
        const int patch_dim_2 = patch_dims[2];
        
        const Real patch_xlo_0 = Real(patch_xlo[0]);
        const Real patch_xlo_1 = Real(patch_xlo[1]);
        const Real patch_xlo_2 = Real(patch_xlo[2]);
        
        const Real dx_0 = Real(dx[0]);
        const Real dx_1 = Real(dx[1]);
        const Real dx_2 = Real(dx[2]);
        
        // Loop over the refine boxes to tag cells for refinement.
        for (int ib = 0; ib < static_cast<int>(d_refine_boxes.size()); ib++)
        {
            // Get the lower and upper coordinates of the refine box.
            const std::vector<Real>& refine_box_x_lo = d_refine_boxes[ib].x_lo;
            const std::vector<Real>& refine_box_x_hi = d_refine_boxes[ib].x_hi;
            
            const Real refine_box_x_lo_0 = refine_box_x_lo[0];
            const Real refine_box_x_lo_1 = refine_box_x_lo[1];
            const Real refine_box_x_lo_2 = refine_box_x_lo[2];
            const Real refine_box_x_hi_0 = refine_box_x_hi[0];
            const Real refine_box_x_hi_1 = refine_box_x_hi[1];
            const Real refine_box_x_hi_2 = refine_box_x_hi[2];
            
            // Get the number of refine levels.
            const int num_refine_levels = d_refine_boxes[ib].num_refine_levels;
            
            if (patch_level < num_refine_levels)
            {
                for (int k = 0; k < patch_dim_2; k++)
                {
                    for (int j = 0; j < patch_dim_1; j++)
                    {
                        HAMERS_PRAGMA_SIMD
                        for (int i = 0; i < patch_dim_0; i++)
                        {
                            // Compute index into linear data array.
                            int idx_cell = i + j*patch_dim_0 + k*patch_dim_0*patch_dim_1;
                            
                            // Compute the coordinates.
                            Real x[3];
                            x[0] = patch_xlo_0 + (Real(i) + half)*dx_0;
                            x[1] = patch_xlo_1 + (Real(j) + half)*dx_1;
                            x[2] = patch_xlo_2 + (Real(k) + half)*dx_2;
                            
                            if (x[0] >= refine_box_x_lo_0 && x[0] <= refine_box_x_hi_0 &&
                                x[1] >= refine_box_x_lo_1 && x[1] <= refine_box_x_hi_1 &&
                                x[2] >= refine_box_x_lo_2 && x[2] <= refine_box_x_hi_2)
                            {
                                tag_ptr[idx_cell] = 1;
                            }
                        }
                    }
                }
            }
        }
    }
}


/*
 * Tag cells on a patch for refinement using refine triangles.
 */
void
RefineRegionsTagger::tagCellsOnPatchUsingRefineTriangles(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags) const
{
    // Get the patch level and path dimensions.
    const int patch_level = patch.getPatchLevelNumber();
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(patch_geom);
#endif
    
    const double* const dx = patch_geom->getDx();
    const double* const patch_xlo = patch_geom->getXLower();
    
    // Get the dimensions of box that covers the interior of Patch.
    hier::Box patch_box = patch.getBox();
    const hier::IntVector patch_dims = patch_box.numberCells();
    
    // Get the pointer to the tags.
    int* tag_ptr = tags->getPointer(0);
    
    const Real half = Real(1)/Real(2);
    
    if (d_dim == tbox::Dimension(1))
    {
        // One-dimensional problem is not supported.
        TBOX_ERROR(d_object_name
            << ": "
            << "Cannot tag cells on patch for one-dimensional problem using refine triangles."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        const int patch_dim_0 = patch_dims[0];
        const int patch_dim_1 = patch_dims[1];
        
        const Real patch_xlo_0 = Real(patch_xlo[0]);
        const Real patch_xlo_1 = Real(patch_xlo[1]);
        
        const Real dx_0 = Real(dx[0]);
        const Real dx_1 = Real(dx[1]);
        
        // Loop over the refine traingles to tag cells for refinement.
        for (int it = 0; it < d_num_refine_triangles; it++)
        {
            // Get the point coordinates of the refine triangle.
            const std::vector<Real>& point_coord_0 = d_refine_triangles[it].point_coord_0;
            const std::vector<Real>& point_coord_1 = d_refine_triangles[it].point_coord_1;
            const std::vector<Real>& point_coord_2 = d_refine_triangles[it].point_coord_2;
            
            const Real x_0 = point_coord_0[0];
            const Real y_0 = point_coord_0[1];
            const Real x_1 = point_coord_1[0];
            const Real y_1 = point_coord_1[1];
            const Real x_2 = point_coord_2[0];
            const Real y_2 = point_coord_2[1];
            
            // Get the number of refine levels.
            const int num_refine_levels = d_refine_triangles[it].num_refine_levels;
            
            // Compute the area of the refine triangle.
            const Real area_refine_triangle = computeTriangleArea(x_0, y_0, x_1, y_1, x_2, y_2);
            
            if (patch_level < num_refine_levels)
            {
                for (int j = 0; j < patch_dim_1; j++)
                {
                    HAMERS_PRAGMA_SIMD
                    for (int i = 0; i < patch_dim_0; i++)
                    {
                        // Compute index into linear data array.
                        int idx_cell = i + j*patch_dim_0;
                        
                        // Compute the coordinates.
                        Real x[2];
                        x[0] = patch_xlo_0 + (Real(i) + half)*dx_0;
                        x[1] = patch_xlo_1 + (Real(j) + half)*dx_1;
                        
                        const Real area_1 = computeTriangleArea(x[0], x[1], x_1, y_1, x_2, y_2);
                        const Real area_2 = computeTriangleArea(x_0, y_0, x[0], x[1], x_2, y_2);
                        const Real area_3 = computeTriangleArea(x_0, y_0, x_1, y_1, x[0], x[1]);
                        
                        // Tag if the point is inside the refine triangle.
                        if (std::abs(area_1 + area_2 + area_3 - area_refine_triangle)/area_refine_triangle
                            < Real(1000)*HAMERS_EPSILON)
                        {
                            tag_ptr[idx_cell] = 1;
                        }
                    }
                }
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        const int patch_dim_0 = patch_dims[0];
        const int patch_dim_1 = patch_dims[1];
        const int patch_dim_2 = patch_dims[2];
        
        const Real patch_xlo_0 = Real(patch_xlo[0]);
        const Real patch_xlo_1 = Real(patch_xlo[1]);
        const Real patch_xlo_2 = Real(patch_xlo[2]);
        
        const Real dx_0 = Real(dx[0]);
        const Real dx_1 = Real(dx[1]);
        const Real dx_2 = Real(dx[2]);
        
        // Loop over the refine traingles to tag cells for refinement.
        for (int it = 0; it < d_num_refine_triangles; it++)
        {
            // Get the point coordinates of the refine triangle.
            const std::vector<Real>& point_coord_0 = d_refine_triangles[it].point_coord_0;
            const std::vector<Real>& point_coord_1 = d_refine_triangles[it].point_coord_1;
            const std::vector<Real>& point_coord_2 = d_refine_triangles[it].point_coord_2;
            
            const Real x_0 = point_coord_0[0];
            const Real y_0 = point_coord_0[1];
            const Real z_0 = point_coord_0[2];
            const Real x_1 = point_coord_1[0];
            const Real y_1 = point_coord_1[1];
            const Real z_1 = point_coord_1[2];
            const Real x_2 = point_coord_2[0];
            const Real y_2 = point_coord_2[1];
            const Real z_2 = point_coord_2[2];
            
            // Get the number of refine levels.
            const int num_refine_levels = d_refine_triangles[it].num_refine_levels;
            
            // Compute the area of the refine triangle.
            Real area_refine_triangle = Real(0);
            
            if (d_refine_triangles[it].direction == DIRECTION::TYPE::X_DIRECTION)
            {
                area_refine_triangle = computeTriangleArea(y_0, z_0, y_1, z_1, y_2, z_2);
            }
            else if (d_refine_triangles[it].direction == DIRECTION::TYPE::Y_DIRECTION)
            {
                area_refine_triangle = computeTriangleArea(x_0, z_0, x_1, z_1, x_2, z_2);
            }
            else if (d_refine_triangles[it].direction == DIRECTION::TYPE::Z_DIRECTION)
            {
                area_refine_triangle = computeTriangleArea(x_0, y_0, x_1, y_1, x_2, y_2);
            }
            
            if (patch_level < num_refine_levels)
            {
                if (d_refine_triangles[it].direction == DIRECTION::TYPE::X_DIRECTION)
                {
                    for (int k = 0; k < patch_dim_2; k++)
                    {
                        for (int j = 0; j < patch_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < patch_dim_0; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dim_0 + k*patch_dim_0*patch_dim_1;
                                
                                // Compute the coordinates.
                                Real x[3];
                                x[0] = patch_xlo_0 + (Real(i) + half)*dx_0;
                                x[1] = patch_xlo_1 + (Real(j) + half)*dx_1;
                                x[2] = patch_xlo_2 + (Real(k) + half)*dx_2;
                                
                                const Real area_1 = computeTriangleArea(x[1], x[2], y_1, z_1, y_2, z_2);
                                const Real area_2 = computeTriangleArea(y_0, z_0, x[1], x[2], y_2, z_2);
                                const Real area_3 = computeTriangleArea(y_0, z_0, y_1, z_1, x[1], x[2]);
                                
                                // Tag if the point is inside the refine triangle.
                                if (std::abs(area_1 + area_2 + area_3 - area_refine_triangle)/area_refine_triangle
                                    < Real(1000)*HAMERS_EPSILON)
                                {
                                    tag_ptr[idx_cell] = 1;
                                }
                            }
                        }
                    }
                }
                else if (d_refine_triangles[it].direction == DIRECTION::TYPE::Y_DIRECTION)
                {
                    for (int k = 0; k < patch_dim_2; k++)
                    {
                        for (int j = 0; j < patch_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < patch_dim_0; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dim_0 + k*patch_dim_0*patch_dim_1;
                                
                                // Compute the coordinates.
                                Real x[3];
                                x[0] = patch_xlo_0 + (Real(i) + half)*dx_0;
                                x[1] = patch_xlo_1 + (Real(j) + half)*dx_1;
                                x[2] = patch_xlo_2 + (Real(k) + half)*dx_2;
                                
                                const Real area_1 = computeTriangleArea(x[0], x[2], x_1, z_1, x_2, z_2);
                                const Real area_2 = computeTriangleArea(x_0, z_0, x[0], x[2], x_2, z_2);
                                const Real area_3 = computeTriangleArea(x_0, z_0, x_1, z_1, x[0], x[2]);
                                
                                // Tag if the point is inside the refine triangle.
                                if (std::abs(area_1 + area_2 + area_3 - area_refine_triangle)/area_refine_triangle
                                    < Real(1000)*HAMERS_EPSILON)
                                {
                                    tag_ptr[idx_cell] = 1;
                                }
                            }
                        }
                    }
                }
                else if (d_refine_triangles[it].direction == DIRECTION::TYPE::Z_DIRECTION)
                {
                    for (int k = 0; k < patch_dim_2; k++)
                    {
                        for (int j = 0; j < patch_dim_1; j++)
                        {
                            HAMERS_PRAGMA_SIMD
                            for (int i = 0; i < patch_dim_0; i++)
                            {
                                // Compute index into linear data array.
                                int idx_cell = i + j*patch_dim_0 + k*patch_dim_0*patch_dim_1;
                                
                                // Compute the coordinates.
                                Real x[3];
                                x[0] = patch_xlo_0 + (Real(i) + half)*dx_0;
                                x[1] = patch_xlo_1 + (Real(j) + half)*dx_1;
                                x[2] = patch_xlo_2 + (Real(k) + half)*dx_2;
                                
                                const Real area_1 = computeTriangleArea(x[0], x[1], x_1, y_1, x_2, y_2);
                                const Real area_2 = computeTriangleArea(x_0, y_0, x[0], x[1], x_2, y_2);
                                const Real area_3 = computeTriangleArea(x_0, y_0, x_1, y_1, x[0], x[1]);
                                
                                // Tag if the point is inside the refine triangle.
                                if (std::abs(area_1 + area_2 + area_3 - area_refine_triangle)/area_refine_triangle
                                    < Real(1000)*HAMERS_EPSILON)
                                {
                                    tag_ptr[idx_cell] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
