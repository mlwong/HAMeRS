#include "extn/variable_fill_patterns/LayerCellVariableFillPattern.hpp"

#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/make_shared.hpp"

const std::string LayerCellVariableFillPattern::s_name_id =
    "LAYER_CELL_FILL_PATTERN";

/*
 *************************************************************************
 *
 * Constructor
 *
 *************************************************************************
 */
LayerCellVariableFillPattern::LayerCellVariableFillPattern(
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
LayerCellVariableFillPattern::~LayerCellVariableFillPattern()
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
LayerCellVariableFillPattern::calculateOverlap(
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
LayerCellVariableFillPattern::getStencilWidth()
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
LayerCellVariableFillPattern::getPatternName() const
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
LayerCellVariableFillPattern::computeStencilBoxes(
    SAMRAI::hier::BoxContainer& stencil_boxes,
    const SAMRAI::hier::Box& dst_box) const
{
    TBOX_ASSERT(stencil_boxes.size() == 0);
    TBOX_ASSERT(dst_box.getDim() == d_dim);
    
    SAMRAI::hier::Box dst_box_expanded(
        SAMRAI::hier::Box::grow(dst_box,
            d_layer_offset));
    
    SAMRAI::hier::Box ghost_box(
        SAMRAI::hier::Box::grow(dst_box_expanded,
            d_layer_width));
    
    stencil_boxes.removeIntersections(ghost_box, dst_box_expanded);
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
LayerCellVariableFillPattern::computeFillBoxesOverlap(
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
