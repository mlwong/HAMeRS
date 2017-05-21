/*************************************************************************
 *
 * This file is modified from FlattenedHierarchy.cpp of the SAMRAI 3.11.2
 * distribution.  For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 ************************************************************************/
#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/HierarchyNeighbors.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"

#include "boost/make_shared.hpp"

/*
 ***************************************************************************
 * Constructor evaluates the hierarchy and fills the visible boxes container.
 ***************************************************************************
 */

ExtendedFlattenedHierarchy::ExtendedFlattenedHierarchy(
    const SAMRAI::hier::PatchHierarchy& hierarchy,
    int coarsest_level,
    int finest_level):
        d_coarsest_level(coarsest_level),
        d_finest_level(finest_level),
        d_patch_hierarchy(&hierarchy)
{
    int num_levels = hierarchy.getNumberOfLevels();
    TBOX_ASSERT(coarsest_level >= 0);
    TBOX_ASSERT(coarsest_level <= finest_level);
    TBOX_ASSERT(finest_level < num_levels);
    
    d_visible_boxes.resize(num_levels);
    d_overlapped_visible_boxes.resize(num_levels);
    
    SAMRAI::hier::LocalId local_id(0);
    SAMRAI::hier::LocalId local_id_overlapped(0);
    
    for (int ln = finest_level; ln >= coarsest_level; ln--)
    {
        const boost::shared_ptr<SAMRAI::hier::PatchLevel>& current_level =
            hierarchy.getPatchLevel(ln);
        
        if (ln != finest_level)
        {
            const SAMRAI::hier::Connector& coarse_to_fine =
                current_level->findConnector(
                    *(hierarchy.getPatchLevel(ln+1)),
                    SAMRAI::hier::IntVector::getOne(hierarchy.getDim()),
                    SAMRAI::hier::CONNECTOR_IMPLICIT_CREATION_RULE,
                    true);
            
            const SAMRAI::hier::Connector& same_level_overlap =
                current_level->findConnector(
                    *(hierarchy.getPatchLevel(ln)),
                    SAMRAI::hier::IntVector::getOne(hierarchy.getDim()),
                    SAMRAI::hier::CONNECTOR_IMPLICIT_CREATION_RULE,
                    true);
            
            const SAMRAI::hier::IntVector& connector_ratio = coarse_to_fine.getRatio();
            
            for (SAMRAI::hier::PatchLevel::iterator ip(current_level->begin());
                 ip != current_level->end();
                 ip++)
            {
                const boost::shared_ptr<SAMRAI::hier::Patch>& patch = *ip;
                const SAMRAI::hier::Box& box = patch->getBox();
                const SAMRAI::hier::BlockId& block_id = box.getBlockId();
                const SAMRAI::hier::BoxId& box_id = box.getBoxId();
                SAMRAI::hier::BoxContainer& visible_boxes = d_visible_boxes[ln][box_id];
                SAMRAI::hier::BoxContainer& overlapped_visible_boxes =
                    d_overlapped_visible_boxes[ln][box_id];
                
                SAMRAI::hier::BoxContainer coarse_boxes(box);
                
                SAMRAI::hier::BoxContainer fine_nbr_boxes;
                if (coarse_to_fine.hasNeighborSet(box_id))
                {
                    coarse_to_fine.getNeighborBoxes(box_id, fine_nbr_boxes);
                }
                if (!fine_nbr_boxes.empty())
                {
                    SAMRAI::hier::BoxContainer fine_boxes;
                    for (SAMRAI::hier::RealBoxConstIterator nbr_itr = fine_nbr_boxes.realBegin();
                         nbr_itr != fine_nbr_boxes.realEnd();
                         nbr_itr++)
                    {
                        if (nbr_itr->getBlockId() == block_id)
                        {
                            fine_boxes.pushBack(*nbr_itr);
                        }
                    }
                    
                    fine_boxes.coarsen(connector_ratio);
                    
                    coarse_boxes.removeIntersections(fine_boxes);
                    coarse_boxes.coalesce();
                }
                
                for (SAMRAI::hier::BoxContainer::iterator itr = coarse_boxes.begin();
                     itr != coarse_boxes.end();
                     itr++)
                {
                    SAMRAI::hier::Box new_box(*itr, local_id, box_id.getOwnerRank());
                    local_id++;
                    visible_boxes.insert(visible_boxes.end(), new_box);
                }
                
                SAMRAI::hier::BoxContainer nbr_boxes;
                if (same_level_overlap.hasNeighborSet(box_id))
                {
                    same_level_overlap.getNeighborBoxes(box_id, nbr_boxes);
                    
                    SAMRAI::hier::BoxContainer overlapped_boxes(coarse_boxes);
                    overlapped_boxes.intersectBoxes(nbr_boxes);
                    
                    for (SAMRAI::hier::BoxContainer::iterator itr = overlapped_boxes.begin();
                         itr != overlapped_boxes.end();
                         itr++)
                    {
                        const SAMRAI::hier::BoxId& nbr_box_id = (*itr).getBoxId();
                        if (nbr_box_id != box_id)
                        {
                            SAMRAI::hier::Box overlapped_box(*itr, local_id_overlapped, box_id.getOwnerRank());
                            local_id_overlapped++;
                            overlapped_visible_boxes.insert(overlapped_visible_boxes.end(), overlapped_box);
                        }
                    }
                }
            }
        }
        else
        {
            const SAMRAI::hier::Connector& same_level_overlap =
                current_level->findConnector(
                    *(hierarchy.getPatchLevel(ln)),
                    SAMRAI::hier::IntVector::getOne(hierarchy.getDim()),
                    SAMRAI::hier::CONNECTOR_IMPLICIT_CREATION_RULE,
                    true);
            
            for (SAMRAI::hier::PatchLevel::iterator ip(current_level->begin());
                 ip != current_level->end();
                 ip++)
            {
                const boost::shared_ptr<SAMRAI::hier::Patch>& patch = *ip;
                const SAMRAI::hier::Box& box = patch->getBox();
                const SAMRAI::hier::BoxId& box_id = box.getBoxId();
                SAMRAI::hier::BoxContainer& visible_boxes = d_visible_boxes[ln][box_id];
                SAMRAI::hier::BoxContainer& overlapped_visible_boxes =
                    d_overlapped_visible_boxes[ln][box_id];
                
                SAMRAI::hier::Box new_box(box, local_id, box.getOwnerRank());
                local_id++;
                visible_boxes.insert(visible_boxes.end(), new_box);
                
                SAMRAI::hier::BoxContainer nbr_boxes;
                if (same_level_overlap.hasNeighborSet(box_id))
                {
                    same_level_overlap.getNeighborBoxes(box_id, nbr_boxes);
                    
                    SAMRAI::hier::BoxContainer overlapped_boxes(new_box);
                    overlapped_boxes.intersectBoxes(nbr_boxes);
                    
                    for (SAMRAI::hier::BoxContainer::iterator itr = overlapped_boxes.begin();
                         itr != overlapped_boxes.end();
                         itr++)
                    {
                        const SAMRAI::hier::BoxId& nbr_box_id = (*itr).getBoxId();
                        if (nbr_box_id != box_id)
                        {
                            SAMRAI::hier::Box overlapped_box(*itr, local_id_overlapped, box_id.getOwnerRank());
                            local_id_overlapped++;
                            overlapped_visible_boxes.insert(overlapped_visible_boxes.end(), overlapped_box);
                        }
                    }
                }
            }
        }
    }
}


/*
 **************************************************************************
 * Destructor
 **************************************************************************
 */

ExtendedFlattenedHierarchy::~ExtendedFlattenedHierarchy()
{
}
