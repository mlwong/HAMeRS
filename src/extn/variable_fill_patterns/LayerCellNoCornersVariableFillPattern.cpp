#include "extn/variable_fill_patterns/LayerCellNoCornersVariableFillPattern.hpp"

#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/make_shared.hpp"

const std::string LayerCellNoCornersVariableFillPattern::s_name_id =
    "LAYER_CELL_NO_CORNERS_FILL_PATTERN";

/*
 *************************************************************************
 *
 * Constructor
 *
 *************************************************************************
 */
LayerCellNoCornersVariableFillPattern::
LayerCellNoCornersVariableFillPattern(
    const SAMRAI::tbox::Dimension& dim,
    const SAMRAI::hier::IntVector& layer_width,
    const SAMRAI::hier::IntVector& layer_offset):
        d_dim(dim),
        d_layer_width(layer_width),
        d_layer_offset(layer_offset)
{
    TBOX_ASSERT(layer_width.getDim() == d_dim);
    TBOX_ASSERT(layer_offset.getDim() == d_dim);
    
    TBOX_ASSERT(layer_width >= SAMRAI::hier::IntVector::getZero(d_dim));
    TBOX_ASSERT(layer_offset >= SAMRAI::hier::IntVector::getZero(d_dim));
}


/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */
LayerCellNoCornersVariableFillPattern::~
LayerCellNoCornersVariableFillPattern()
{
}


/*
 *************************************************************************
 *
 * Calculate the overlap according to the desired pattern
 *
 *************************************************************************
 */
boost::shared_ptr<SAMRAI::hier::BoxOverlap>
LayerCellNoCornersVariableFillPattern::calculateOverlap(
    const SAMRAI::hier::BoxGeometry& dst_geometry,
    const SAMRAI::hier::BoxGeometry& src_geometry,
    const SAMRAI::hier::Box& dst_patch_box,
    const SAMRAI::hier::Box& src_mask,
    const SAMRAI::hier::Box& fill_box,
    const bool overwrite_interior,
    const SAMRAI::hier::Transformation& transformation) const
{
    TBOX_ASSERT_OBJDIM_EQUALITY2(dst_patch_box, src_mask);
    
    SAMRAI::hier::BoxContainer stencil_boxes;
    computeStencilBoxes(stencil_boxes, dst_patch_box);
    
    return dst_geometry.calculateOverlap(
        src_geometry,
        src_mask,
        fill_box,
        overwrite_interior,
        transformation,
        stencil_boxes);
}


/*
 *************************************************************************
 *
 * Return the stencil width
 *
 *************************************************************************
 */
const SAMRAI::hier::IntVector&
LayerCellNoCornersVariableFillPattern::getStencilWidth()
{
    return d_layer_width;
}


/*
 *************************************************************************
 *
 * Return the string name identifier
 *
 *************************************************************************
 */
const std::string&
LayerCellNoCornersVariableFillPattern::getPatternName() const
{
    return s_name_id;
}


/*
 *************************************************************************
 *
 * Compute the boxes for the stencil around a given patch box
 *
 *************************************************************************
 */
void
LayerCellNoCornersVariableFillPattern::computeStencilBoxes(
    SAMRAI::hier::BoxContainer& stencil_boxes,
    const SAMRAI::hier::Box& dst_box) const
{
    TBOX_ASSERT(stencil_boxes.size() == 0);
    TBOX_ASSERT(dst_box.getDim() == d_dim);
    
    for (unsigned short i = 0; i < d_dim.getValue(); i++)
    {
        if (d_layer_width[i] > 0)
        {
            SAMRAI::hier::Box low_box(dst_box);
            low_box.setUpper(i, dst_box.lower(i) - d_layer_offset[i] - 1);
            low_box.setLower(i, low_box.upper(i) - d_layer_width[i] + 1);
            stencil_boxes.pushFront(low_box);
            
            SAMRAI::hier::Box high_box(dst_box);
            high_box.setLower(i, dst_box.upper(i) + d_layer_offset[i] + 1);
            high_box.setUpper(i, high_box.lower(i) + d_layer_width[i] - 1);
            stencil_boxes.pushFront(high_box);
        }
    }
}


/*
 *************************************************************************
 *
 * Compute BoxOverlap that specifies data to be filled by refinement
 * operator.
 *
 *************************************************************************
 */
boost::shared_ptr<SAMRAI::hier::BoxOverlap>
LayerCellNoCornersVariableFillPattern::computeFillBoxesOverlap(
    const SAMRAI::hier::BoxContainer& fill_boxes,
    const SAMRAI::hier::BoxContainer& node_fill_boxes,
    const SAMRAI::hier::Box& patch_box,
    const SAMRAI::hier::Box& data_box,
    const SAMRAI::hier::PatchDataFactory& pdf) const
{
    NULL_USE(pdf);
    NULL_USE(node_fill_boxes);
    
    SAMRAI::hier::BoxContainer stencil_boxes;
    computeStencilBoxes(stencil_boxes, patch_box);
    
    SAMRAI::hier::BoxContainer overlap_boxes(fill_boxes);
    overlap_boxes.intersectBoxes(data_box);
    overlap_boxes.intersectBoxes(stencil_boxes);
    
    return boost::make_shared<SAMRAI::pdat::CellOverlap>(
        overlap_boxes,
        SAMRAI::hier::Transformation(SAMRAI::hier::IntVector::getZero(patch_box.getDim())));
}
