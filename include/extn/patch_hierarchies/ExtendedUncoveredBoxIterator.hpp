/*************************************************************************
 *
 * This file is modified from UncoveredBoxIterator.h of the SAMRAI 3.11.2
 * distribution.  For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 ************************************************************************/

#ifndef EXTENDED_UNCOVERED_BOX_ITERATOR
#define EXTENDED_UNCOVERED_BOX_ITERATOR

#include "HAMeRS_config.hpp"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/PatchHierarchy.h"

#include "boost/shared_ptr.hpp"
#include <utility>
#include <vector>

class ExtendedFlattenedHierarchy;

/*!
 * @brief An iterator over the uncovered Boxes in a hierarhcy.
 *
 * This iterator points to a pair consisting of the current uncovered Box and the Patch with which it is associated.  Note that in the case of overlapping
 * Patches in which the overlap is uncovered, the overlapping region will
 * appear multiple times in the iteration.  The number of appearances is equal
 * to the number of Patches that overlap that region.
 */
class ExtendedUncoveredBoxIterator
{
    friend class SAMRAI::hier::PatchHierarchy;
    friend class ExtendedFlattenedHierarchy;

public:
    /*!
     * @brief Copy constructor.
     *
     * @param[in] other The iterator being copied.
     */
    ExtendedUncoveredBoxIterator(
        const ExtendedUncoveredBoxIterator& other);
    
    /*!
     * Destructor.
     */
    ~ExtendedUncoveredBoxIterator();
    
    /*!
     * @brief Assignment operator.
     *
     * @param[in] rhs The right hand side of the assignment.
     */
    ExtendedUncoveredBoxIterator&
    operator = (
        const ExtendedUncoveredBoxIterator& rhs);
    
    /*!
     * @brief Dereference operator mimicking a pointer dereference.
     */
    const std::pair<boost::shared_ptr<SAMRAI::hier::Patch>, SAMRAI::hier::Box>&
    operator * () const;
    
    /*!
     * @brief Dereference operator mimicking a pointer dereference.
     */
    const std::pair<boost::shared_ptr<SAMRAI::hier::Patch>, SAMRAI::hier::Box> *
    operator -> () const;
    
    /*!
     * @brief Equality comparison.
     *
     * @param[in] rhs The right hand side of the comparison.
     */
    bool
    operator == (
        const ExtendedUncoveredBoxIterator& rhs) const;
    
    /*!
     * @brief Inequality comparison.
     *
     * @param[in] rhs The right hand side of the comparison.
     */
    bool
    operator != (
        const ExtendedUncoveredBoxIterator& rhs) const;
    
    /*!
     * @brief Pre-increment iterator.
     *
     * Pre-increment increment the iterator and returns the incremented
     * state.
     */
    ExtendedUncoveredBoxIterator&
    operator ++ ();
    
    /*!
     * @brief Post-increment iterator.
     *
     * Post-increment saves the iterator, increments it and returns the
     * saved iterator.
     */
    ExtendedUncoveredBoxIterator
    operator ++ (
        int);

private:
    /*!
     * @brief Unimplemented default constructor.
     */
    ExtendedUncoveredBoxIterator();
    
    /*!
     * @brief Constructor.
     *
     * @param[in] hierarchy The hierarchy over which the iteration will occur.
     * @param[in] begin If true the iteration starts from the beginning.
     */
    ExtendedUncoveredBoxIterator(
        const SAMRAI::hier::PatchHierarchy* hierarchy,
        bool begin);
    
    ExtendedUncoveredBoxIterator(
        const ExtendedFlattenedHierarchy* hierarchy,
        bool begin);
    
    /*!
     * @brief Private method to do work common to both pre and post increments.
     */
    void
    incrementIterator();
    
    /*!
     * @brief Private method to find first uncovered box
     */
    void
    findFirstUncoveredBox();
    
    /*!
     * @brief Set the item based on the current patch id and uncovered box.
     */
    void
    setIteratorItem();
    
    /* The PatchHierarchy on which this iterator operates. */
    const SAMRAI::hier::PatchHierarchy* d_hierarchy;
    const ExtendedFlattenedHierarchy* d_extended_flattened_hierarchy;
    
    bool d_allocated_extended_flattened_hierarchy;
    
    /* The current level in the PatchHierarchy. */
    int d_level_num;
    
    /* The id of the current patch */
    SAMRAI::hier::BoxId d_current_patch_id;
    
    /* The iterator over the uncovered boxes for the current patch. */
    SAMRAI::hier::BoxContainer::const_iterator d_uncovered_boxes_itr;
    
    /* The iterator at the end of the uncovered boxes for the current patch. */
    SAMRAI::hier::BoxContainer::const_iterator d_uncovered_boxes_itr_end;
    
    /* The current item in the iteration. */
    std::pair<boost::shared_ptr<SAMRAI::hier::Patch>, SAMRAI::hier::Box>* d_item;
    
    /* The number of the finest level in the hierarchy. */
    int d_finest_level_num;

};

#endif /* EXTENDED_UNCOVERED_BOX_ITERATOR */
