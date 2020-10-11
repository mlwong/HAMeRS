/*************************************************************************
 *
 * This file is modified from UncoveredBoxIterator.c of the SAMRAI 3.11.2
 * distribution.  For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 ************************************************************************/
#include "extn/patch_hierarchies/ExtendedUncoveredBoxIterator.hpp"

#include "extn/patch_hierarchies/ExtendedFlattenedHierarchy.hpp"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"

ExtendedUncoveredBoxIterator::ExtendedUncoveredBoxIterator(
    const SAMRAI::hier::PatchHierarchy* hierarchy,
    bool begin):
    d_hierarchy(hierarchy),
    d_uncovered_boxes_itr(SAMRAI::hier::BoxContainer().begin()),
    d_uncovered_boxes_itr_end(SAMRAI::hier::BoxContainer().end()),
    d_item(0)
{
    TBOX_ASSERT(hierarchy);
 
    d_finest_level_num = d_hierarchy->getFinestLevelNumber();
    if (begin)
    {
        d_level_num = -1;
        d_extended_flattened_hierarchy = new ExtendedFlattenedHierarchy(*d_hierarchy, 0, d_finest_level_num);
        d_allocated_extended_flattened_hierarchy = true;
        findFirstUncoveredBox();
    }
    else
    {
        d_level_num = d_finest_level_num + 1;
        d_extended_flattened_hierarchy = 0;
        d_allocated_extended_flattened_hierarchy = false;
    }
}


ExtendedUncoveredBoxIterator::ExtendedUncoveredBoxIterator(
    const ExtendedFlattenedHierarchy* extended_flattened_hierarchy,
    bool begin):
    d_hierarchy(&(extended_flattened_hierarchy->getPatchHierarchy())),
    d_uncovered_boxes_itr(SAMRAI::hier::BoxContainer().begin()),
    d_uncovered_boxes_itr_end(SAMRAI::hier::BoxContainer().end()),
    d_item(0)
{
    TBOX_ASSERT(extended_flattened_hierarchy);
 
    d_finest_level_num = d_hierarchy->getFinestLevelNumber();
    if (begin)
    {
        d_level_num = -1;
        d_extended_flattened_hierarchy = extended_flattened_hierarchy;
        d_allocated_extended_flattened_hierarchy = false;
        findFirstUncoveredBox();
    }
    else
    {
        d_level_num = d_finest_level_num + 1;
        d_extended_flattened_hierarchy = 0;
        d_allocated_extended_flattened_hierarchy = false;
    }
}


ExtendedUncoveredBoxIterator::ExtendedUncoveredBoxIterator(
    const ExtendedUncoveredBoxIterator& other):
    d_hierarchy(other.d_hierarchy),
    d_extended_flattened_hierarchy(0),
    d_allocated_extended_flattened_hierarchy(false),
    d_level_num(other.d_level_num),
    d_current_patch_id(other.d_current_patch_id),
    d_uncovered_boxes_itr(other.d_uncovered_boxes_itr),
    d_uncovered_boxes_itr_end(other.d_uncovered_boxes_itr_end),
    d_item(0),
    d_finest_level_num(other.d_finest_level_num)
{
    if (other.d_item)
    {
        d_item = new std::pair<HAMERS_SHARED_PTR<SAMRAI::hier::Patch>, SAMRAI::hier::Box>(*other.d_item);
    }
    if (other.d_extended_flattened_hierarchy)
    {
        d_extended_flattened_hierarchy =
            new ExtendedFlattenedHierarchy(*other.d_extended_flattened_hierarchy);
        d_allocated_extended_flattened_hierarchy = true;
        if (d_level_num <= d_finest_level_num)
        {
            TBOX_ASSERT(d_item);
            const SAMRAI::hier::Box& patch_box = d_item->first->getBox();
            const SAMRAI::hier::BoxContainer& visible_boxes =
               d_extended_flattened_hierarchy->getVisibleBoxes(patch_box, d_level_num);
            
            SAMRAI::hier::BoxContainer::const_iterator itr = visible_boxes.begin(); 
               for( ; itr != visible_boxes.end(); ++itr)
            {
                if (itr->getBoxId() == d_item->second.getBoxId() &&
                    itr->isSpatiallyEqual(d_item->second))
                {
                    d_uncovered_boxes_itr = itr;
                    d_uncovered_boxes_itr_end = visible_boxes.end();
                    break;
                }
            }
            
            if (itr == visible_boxes.end())
            {
                d_uncovered_boxes_itr = itr;
                d_uncovered_boxes_itr_end = itr;
            }
        }
    }
}


ExtendedUncoveredBoxIterator::~ExtendedUncoveredBoxIterator()
{
    if (d_item)
    {
        delete d_item;
    }
    if (d_extended_flattened_hierarchy && d_allocated_extended_flattened_hierarchy)
    {
        delete d_extended_flattened_hierarchy;
    } 
}


ExtendedUncoveredBoxIterator&
ExtendedUncoveredBoxIterator::operator = (
    const ExtendedUncoveredBoxIterator& rhs)
{
    if (this != &rhs)
    {
        d_hierarchy = rhs.d_hierarchy;
        d_level_num = rhs.d_level_num;
        d_uncovered_boxes_itr = rhs.d_uncovered_boxes_itr;
        d_uncovered_boxes_itr_end = rhs.d_uncovered_boxes_itr_end;
        d_current_patch_id = rhs.d_current_patch_id; 
        if (d_item)
        {
            delete d_item;
        }
        if (rhs.d_item)
        {
            d_item = new std::pair<HAMERS_SHARED_PTR<SAMRAI::hier::Patch>, SAMRAI::hier::Box>(*rhs.d_item);
            d_item->first = rhs.d_item->first;
            d_item->second = rhs.d_item->second;
        }
        else
        {
            d_item = 0;
        } 
        d_finest_level_num = rhs.d_finest_level_num;
        if (d_extended_flattened_hierarchy && d_allocated_extended_flattened_hierarchy)
        {
            delete d_extended_flattened_hierarchy;
        }
        d_extended_flattened_hierarchy = 0;
        d_allocated_extended_flattened_hierarchy = false;
        if (rhs.d_extended_flattened_hierarchy)
        {
            d_extended_flattened_hierarchy =
                new ExtendedFlattenedHierarchy(*rhs.d_extended_flattened_hierarchy);
            d_allocated_extended_flattened_hierarchy = true;
            if (d_level_num <= d_finest_level_num)
            {
                TBOX_ASSERT(d_item);
                const SAMRAI::hier::Box& patch_box = d_item->first->getBox();
                const SAMRAI::hier::BoxContainer& visible_boxes =
                d_extended_flattened_hierarchy->getVisibleBoxes(patch_box, d_level_num);
                
                SAMRAI::hier::BoxContainer::const_iterator itr = visible_boxes.begin();
                for( ; itr != visible_boxes.end(); ++itr)
                {
                    if (itr->getBoxId() == d_item->second.getBoxId() &&
                        itr->isSpatiallyEqual(d_item->second))
                    {
                        
                        d_uncovered_boxes_itr = itr;
                        d_uncovered_boxes_itr_end = visible_boxes.end();
                        break;
                    }
                }
                
                if (itr == visible_boxes.end())
                {
                    d_uncovered_boxes_itr = itr;
                    d_uncovered_boxes_itr_end = itr;
                }
            }
        }
    }
 
    return *this;
}


const std::pair<HAMERS_SHARED_PTR<SAMRAI::hier::Patch>, SAMRAI::hier::Box>&
ExtendedUncoveredBoxIterator::operator * () const
{
    return *d_item;
}


const std::pair<HAMERS_SHARED_PTR<SAMRAI::hier::Patch>, SAMRAI::hier::Box> *
ExtendedUncoveredBoxIterator::operator -> () const
{
    return d_item;
}


bool
ExtendedUncoveredBoxIterator::operator == (
   const ExtendedUncoveredBoxIterator& rhs) const
{
    // Frist check and see if the iterators are working on the same hierarchies
    // and levels.  If not then they are not equal.
    bool result = d_hierarchy == rhs.d_hierarchy &&
        d_level_num == rhs.d_level_num;
    
    if (d_extended_flattened_hierarchy == 0 && rhs.d_extended_flattened_hierarchy != 0)
    {
        result = false;
    }
    if (d_extended_flattened_hierarchy != 0 && rhs.d_extended_flattened_hierarchy == 0)
    {
        result = false;
    }
    if (d_item == 0 && rhs.d_item != 0)
    {
        result = false;
    }
    if (d_item != 0 && rhs.d_item == 0)
    {
        result = false;
    }
    if (result)
    {
        if (d_item == 0 && rhs.d_item == 0)
        {
            result = true;
        }
        if (d_item && rhs.d_item)
        {
            if (d_item->second.isIdEqual(rhs.d_item->second) &&
                d_item->second.isSpatiallyEqual(rhs.d_item->second) &&
                d_item->first->getBox().isIdEqual(rhs.d_item->first->getBox()) &&
                d_item->first->getBox().isSpatiallyEqual(rhs.d_item->first->getBox()))
            {
                result = true;
            }
            else
            {
                result = false;
            }
        }
    }
 
    return result;
}


bool
ExtendedUncoveredBoxIterator::operator != (
   const ExtendedUncoveredBoxIterator& rhs) const
{
    return !(*this == rhs);
}


ExtendedUncoveredBoxIterator&
ExtendedUncoveredBoxIterator::operator ++ ()
{
    incrementIterator();
    return *this;
}


ExtendedUncoveredBoxIterator
ExtendedUncoveredBoxIterator::operator ++ (int)
{
    // Save the state of the iterator.
    ExtendedUncoveredBoxIterator tmp(*this);
 
    incrementIterator();
 
    // Return iterator in original state.
    return tmp;
}

void
ExtendedUncoveredBoxIterator::incrementIterator()
{
    /*
     * Increment the iterator over the uncovered boxes for the current patch.
     * If we reach the end of the uncovered boxes for the current patch,
     * look for the next patch with uncovered boxes, moving to finer levels if
     * necessary
     */
    ++d_uncovered_boxes_itr;
    if (d_uncovered_boxes_itr != d_uncovered_boxes_itr_end)
    {
        d_current_patch_id = d_item->first->getBox().getBoxId();
        setIteratorItem();
    }
    else
    {
        bool id_found = false;
        
        bool new_level = false;
        while (d_level_num <= d_finest_level_num)
        {
            HAMERS_SHARED_PTR<SAMRAI::hier::PatchLevel> this_level =
                d_hierarchy->getPatchLevel(d_level_num);
            
            const SAMRAI::hier::BoxContainer& this_level_boxes =
                this_level->getBoxLevel()->getBoxes();
            
            for (SAMRAI::hier::RealBoxConstIterator this_itr = this_level_boxes.realBegin();
                 this_itr != this_level_boxes.realEnd(); ++this_itr)
            {
                if (!new_level &&
                    this_itr->getBoxId() <= d_item->first->getBox().getBoxId())
                {
                    continue;
                }
                const SAMRAI::hier::BoxContainer& uncovered_boxes =
                    d_extended_flattened_hierarchy->getVisibleBoxes(*this_itr, d_level_num);
                
                if (!uncovered_boxes.empty())
                {
                    d_uncovered_boxes_itr = uncovered_boxes.begin();
                    d_uncovered_boxes_itr_end = uncovered_boxes.end();
                    d_current_patch_id = this_itr->getBoxId();
                    id_found = true;
                    break;
                }
            }
            if (id_found)
            {
                break;
            }
            else
            {
                ++d_level_num;
                new_level = true;
            }
        }
        
        if (id_found)
        {
            setIteratorItem();
        }
        else
        {
            if (d_extended_flattened_hierarchy && d_allocated_extended_flattened_hierarchy)
            {
                delete d_extended_flattened_hierarchy;
            }
            d_extended_flattened_hierarchy = 0;
            if (d_item)
            {
                delete d_item;
                d_item = 0; 
            }
        }
    }
}


void
ExtendedUncoveredBoxIterator::findFirstUncoveredBox()
{
    ++d_level_num;
    HAMERS_SHARED_PTR<SAMRAI::hier::PatchLevel> this_level =
       d_hierarchy->getPatchLevel(d_level_num);
    
    bool id_found = false;
    
    while (d_level_num <= d_finest_level_num)
    { 
        HAMERS_SHARED_PTR<SAMRAI::hier::PatchLevel> this_level =
            d_hierarchy->getPatchLevel(d_level_num);
        
        const SAMRAI::hier::BoxContainer& this_level_boxes =
            this_level->getBoxLevel()->getBoxes();
        
        for (SAMRAI::hier::RealBoxConstIterator this_itr = this_level_boxes.realBegin();
             this_itr != this_level_boxes.realEnd(); ++this_itr)
        {
            const SAMRAI::hier::BoxContainer& uncovered_boxes =
                d_extended_flattened_hierarchy->getVisibleBoxes(*this_itr, d_level_num);
            
            if (!uncovered_boxes.empty())
            {
                d_uncovered_boxes_itr = uncovered_boxes.begin();
                d_uncovered_boxes_itr_end = uncovered_boxes.end();
                d_current_patch_id = this_itr->getBoxId(); 
                id_found = true;
                break;
            }
        }
        if (id_found)
        {
            break;
        }
        else
        {
            ++d_level_num;
        }
    }
 
    if (id_found)
    {
        setIteratorItem();
    }
    else
    {
        if (d_item)
        {
            delete d_item;
            d_item = 0;
        }
        if (d_extended_flattened_hierarchy && d_allocated_extended_flattened_hierarchy)
        {
            delete d_extended_flattened_hierarchy;
        }
        d_extended_flattened_hierarchy = 0;
        if (d_item)
        {
            delete d_item;
            d_item = 0;
        }
    }
}


void
ExtendedUncoveredBoxIterator::setIteratorItem()
{
    // Get the current uncovered box.
    const SAMRAI::hier::Box& cur_box = *d_uncovered_boxes_itr;
    HAMERS_SHARED_PTR<SAMRAI::hier::PatchLevel> this_level =
        d_hierarchy->getPatchLevel(d_level_num);
 
    // Update item with the current originating patch and the current box.
    if (d_item)
    {
        d_item->first =
            this_level->getPatch(d_current_patch_id);
        d_item->second = cur_box;
    }
    else
    {
        d_item =
            new std::pair<HAMERS_SHARED_PTR<SAMRAI::hier::Patch>, SAMRAI::hier::Box>(
                this_level->getPatch(d_current_patch_id),
                cur_box);
    }
}
