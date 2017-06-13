/*************************************************************************
 *
 * This file is modified from RefineAlgorithm.C of the SAMRAI version 3.11.2
 * distribution.  For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of the SAMRAI distribution.
 *
 ************************************************************************/
#include "extn/transfer/ExtendedRefineAlgorithm.hpp"

#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/xfer/StandardRefineTransactionFactory.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/make_shared.hpp"

/*
 **************************************************************************************************
 *
 * Default constructor creates a new RefineClasses object.
 *
 **************************************************************************************************
 */
ExtendedRefineAlgorithm::ExtendedRefineAlgorithm():
    d_refine_classes(boost::make_shared<SAMRAI::xfer::RefineClasses>()),
    d_schedule_created(false)
{
}


/*
 **************************************************************************************************
 *
 * The destructor implicitly deletes the list storage associated with the refine algorithm.
 *
 **************************************************************************************************
 */
ExtendedRefineAlgorithm::~ExtendedRefineAlgorithm()
{
}


/*
 **************************************************************************************************
 *
 * Register a refine operation that will not require time interpolation.
 *
 **************************************************************************************************
 */
void
ExtendedRefineAlgorithm::registerRefine(
    const int dst,
    const int src,
    const int scratch,
    const boost::shared_ptr<SAMRAI::hier::RefineOperator>& oprefine,
    const boost::shared_ptr<SAMRAI::xfer::VariableFillPattern>& var_fill_pattern)
{
    if (d_schedule_created)
    {
        TBOX_ERROR("ExtendedRefineAlgorithm::registerRefine error..."
            << "\nCannot call registerRefine with a ExtendedRefineAlgorithm"
            << "\nobject that has already been used to create a schedule."
            << std::endl);
    }
    
    SAMRAI::xfer::RefineClasses::Data data;
    
    data.d_dst = dst;
    data.d_src = src;
    data.d_src_told = -1;
    data.d_src_tnew = -1;
    data.d_scratch = scratch;
    data.d_fine_bdry_reps_var = SAMRAI::hier::VariableDatabase::getDatabase()->
       getPatchDescriptor()->getPatchDataFactory(dst)->
       fineBoundaryRepresentsVariable();
    data.d_time_interpolate = false;
    data.d_oprefine = oprefine;
    data.d_optime.reset();
    data.d_tag = -1;
    if (var_fill_pattern)
    {
        data.d_var_fill_pattern = var_fill_pattern;
    }
    else
    {
        data.d_var_fill_pattern.reset(new SAMRAI::xfer::BoxGeometryVariableFillPattern());
    }
    
    d_refine_classes->insertEquivalenceClassItem(data);
}


/*
 **************************************************************************************************
 *
 * Register a refine operation that will require time interpolation.
 *
 **************************************************************************************************
 */
void
ExtendedRefineAlgorithm::registerRefine(
    const int dst,
    const int src,
    const int src_told,
    const int src_tnew,
    const int scratch,
    const boost::shared_ptr<SAMRAI::hier::RefineOperator>& oprefine,
    const boost::shared_ptr<SAMRAI::hier::TimeInterpolateOperator>& optime,
    const boost::shared_ptr<SAMRAI::xfer::VariableFillPattern>& var_fill_pattern)
{
    TBOX_ASSERT(optime);
    
    if (d_schedule_created)
    {
        TBOX_ERROR("ExtendedRefineAlgorithm::registerRefine error..."
            << "\nCannot call registerRefine with a ExtendedRefineAlgorithm object"
            << "\nthat has already been used to create a schedule."
            << std::endl);
    }
    
    SAMRAI::xfer::RefineClasses::Data data;
    
    data.d_dst = dst;
    data.d_src = src;
    data.d_src_told = src_told;
    data.d_src_tnew = src_tnew;
    data.d_scratch = scratch;
    data.d_fine_bdry_reps_var = SAMRAI::hier::VariableDatabase::getDatabase()->
       getPatchDescriptor()->getPatchDataFactory(dst)->
       fineBoundaryRepresentsVariable();
    data.d_time_interpolate = true;
    data.d_oprefine = oprefine;
    data.d_optime = optime;
    data.d_tag = -1;
    if (var_fill_pattern)
    {
        data.d_var_fill_pattern = var_fill_pattern;
    }
    else
    {
        data.d_var_fill_pattern.reset(new SAMRAI::xfer::BoxGeometryVariableFillPattern());
    }
    
    d_refine_classes->insertEquivalenceClassItem(data);
}


/*
 **************************************************************************************************
 *
 * Create a communication schedule that will copy data from the interiors of the specified level
 * into the ghost cells and interiors of the same level.
 *
 **************************************************************************************************
 */
boost::shared_ptr<ExtendedRefineSchedule>
ExtendedRefineAlgorithm::createSchedule(
   const boost::shared_ptr<SAMRAI::hier::PatchLevel>& level,
   SAMRAI::xfer::RefinePatchStrategy* patch_strategy,
   const boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory>& transaction_factory)
{
    TBOX_ASSERT(level);
    
    d_schedule_created = true;
    
    boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory> trans_factory(
        transaction_factory);
    
    if (!trans_factory)
    {
        trans_factory.reset(new SAMRAI::xfer::StandardRefineTransactionFactory);
    }
    
    boost::shared_ptr<SAMRAI::xfer::PatchLevelFullFillPattern> fill_pattern(
        boost::make_shared<SAMRAI::xfer::PatchLevelFullFillPattern>());
    
    return boost::make_shared<ExtendedRefineSchedule>(
        fill_pattern,
        level,
        level,
        d_refine_classes,
        trans_factory,
        patch_strategy);
}


/*
 **************************************************************************************************
 *
 * Create a communication schedule that will copy data from the interiors of the specified level
 * into the ghost cells and interiors of the same level.
 *
 **************************************************************************************************
 */
boost::shared_ptr<ExtendedRefineSchedule>
ExtendedRefineAlgorithm::createSchedule(
   const boost::shared_ptr<SAMRAI::xfer::PatchLevelFillPattern>& fill_pattern,
   const boost::shared_ptr<SAMRAI::hier::PatchLevel>& level,
   SAMRAI::xfer::RefinePatchStrategy* patch_strategy,
   const boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory>& transaction_factory)
{
    TBOX_ASSERT(level);
    
    d_schedule_created = true;
    
    boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory> trans_factory(
       transaction_factory);
    
    if (!trans_factory)
    {
       trans_factory.reset(new SAMRAI::xfer::StandardRefineTransactionFactory);
    }
    
    return boost::make_shared<ExtendedRefineSchedule>(
        fill_pattern,
        level,
        level,
        d_refine_classes,
        trans_factory,
        patch_strategy);
}


/*
 **************************************************************************************************
 *
 * Create a communication schedule that will copy data from the interiors of the source level into
 * the ghost cell and interiors of the destination level.
 *
 **************************************************************************************************
 */
boost::shared_ptr<ExtendedRefineSchedule>
ExtendedRefineAlgorithm::createSchedule(
   const boost::shared_ptr<SAMRAI::hier::PatchLevel>& dst_level,
   const boost::shared_ptr<SAMRAI::hier::PatchLevel>& src_level,
   SAMRAI::xfer::RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory>& transaction_factory)
{
    // TBOX_ERROR("Untried method!  I think this method should work, but it's never been excercised.
    // When code crashes here, remove this line and rerun.  If problem continues, it could well be
    // due to excercising this code.  --BTNG");
    
    TBOX_ASSERT(dst_level);
    TBOX_ASSERT(src_level);
    TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *src_level);
    
    d_schedule_created = true;
    
    boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory> trans_factory(
        transaction_factory);
    
    if (!trans_factory)
    {
        trans_factory.reset(new SAMRAI::xfer::StandardRefineTransactionFactory);
    }
    
    boost::shared_ptr<SAMRAI::xfer::PatchLevelFullFillPattern> fill_pattern(
        boost::make_shared<SAMRAI::xfer::PatchLevelFullFillPattern>());
    
    return boost::make_shared<ExtendedRefineSchedule>(
        fill_pattern,
        dst_level,
        src_level,
        d_refine_classes,
        trans_factory,
        patch_strategy,
        use_time_refinement);
}


/*
 **************************************************************************************************
 *
 * Create a communication schedule that will copy data from the interiors of the source level into
 * the ghost cell and interiors of the destination level.
 *
 **************************************************************************************************
 */
boost::shared_ptr<ExtendedRefineSchedule>
ExtendedRefineAlgorithm::createSchedule(
    const boost::shared_ptr<SAMRAI::xfer::PatchLevelFillPattern>& fill_pattern,
    const boost::shared_ptr<SAMRAI::hier::PatchLevel>& dst_level,
    const boost::shared_ptr<SAMRAI::hier::PatchLevel>& src_level,
    SAMRAI::xfer::RefinePatchStrategy* patch_strategy,
    bool use_time_refinement,
    const boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory>& transaction_factory)
{
    TBOX_ASSERT(dst_level);
    TBOX_ASSERT(src_level);
    TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *src_level);
    
    d_schedule_created = true;
    
    boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory> trans_factory(
        transaction_factory);
    
    if (!trans_factory)
    {
        trans_factory.reset(new SAMRAI::xfer::StandardRefineTransactionFactory);
    }
    
    return boost::make_shared<ExtendedRefineSchedule>(
        fill_pattern,
        dst_level,
        src_level,
        d_refine_classes,
        trans_factory,
        patch_strategy,
        use_time_refinement);
}


/*
 **************************************************************************************************
 *
 * Create a communication schedule that copies data from the interiors of the same level and coarser
 * levels into the interior and boundary cells of the given level.
 *
 **************************************************************************************************
 */
boost::shared_ptr<ExtendedRefineSchedule>
ExtendedRefineAlgorithm::createSchedule(
    const boost::shared_ptr<SAMRAI::hier::PatchLevel>& level,
    const int next_coarser_level,
    const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    SAMRAI::xfer::RefinePatchStrategy* patch_strategy,
    bool use_time_refinement,
    const boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory>& transaction_factory)
{
    // Do we all agree on the destination box_level?
    TBOX_ASSERT(level);
    TBOX_ASSERT((next_coarser_level == -1) || hierarchy);
    
#ifdef HAMERS_DEBUG_CHECK_DIM_ASSERTIONS
    if (hierarchy)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY2(*level, *hierarchy);
    }
#endif
    
    d_schedule_created = true;
    
    boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory> trans_factory(
       transaction_factory);
    
    if (!trans_factory)
    {
        trans_factory.reset(new SAMRAI::xfer::StandardRefineTransactionFactory);
    }
    
    boost::shared_ptr<SAMRAI::xfer::PatchLevelFullFillPattern> fill_pattern(
        boost::make_shared<SAMRAI::xfer::PatchLevelFullFillPattern>());
    
    return boost::make_shared<ExtendedRefineSchedule>(
        fill_pattern,
        level,
        level,
        next_coarser_level,
        hierarchy,
        d_refine_classes,
        trans_factory,
        patch_strategy,
        use_time_refinement);
}


/*
 **************************************************************************************************
 *
 * Create a communication schedule that copies data from the interiors of the same level and coarser
 * levels into the interior and boundary cells of the given level.
 *
 **************************************************************************************************
 */
boost::shared_ptr<ExtendedRefineSchedule>
ExtendedRefineAlgorithm::createSchedule(
    const boost::shared_ptr<SAMRAI::xfer::PatchLevelFillPattern>& fill_pattern,
    const boost::shared_ptr<SAMRAI::hier::PatchLevel>& level,
    const int next_coarser_level,
    const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    SAMRAI::xfer::RefinePatchStrategy* patch_strategy,
    bool use_time_refinement,
    const boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory>& transaction_factory)
{
   // Do we all agree on the destination box_level?
   TBOX_ASSERT(level);
   TBOX_ASSERT((next_coarser_level == -1) || hierarchy);
   
#ifdef HAMERS_DEBUG_CHECK_DIM_ASSERTIONS
    if (hierarchy)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY2(*level, *hierarchy);
    }
#endif
    
    d_schedule_created = true;
    
    boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory> trans_factory(
       transaction_factory);
    
    if (!trans_factory)
    {
        trans_factory.reset(new SAMRAI::xfer::StandardRefineTransactionFactory);
    }
    
    return boost::make_shared<ExtendedRefineSchedule>(
        fill_pattern,
        level,
        level,
        next_coarser_level,
        hierarchy,
        d_refine_classes,
        trans_factory,
        patch_strategy,
        use_time_refinement);
}


/*
 **************************************************************************************************
 *
 * Create a communication schedule that copies data from the interiors of the old level and coarser
 * levels into the ghost cells and interior cells of the given new level.
 *
 **************************************************************************************************
 */
boost::shared_ptr<ExtendedRefineSchedule>
ExtendedRefineAlgorithm::createSchedule(
   const boost::shared_ptr<SAMRAI::hier::PatchLevel>& dst_level,
   const boost::shared_ptr<SAMRAI::hier::PatchLevel>& src_level,
   const int next_coarser_level,
   const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
   SAMRAI::xfer::RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory>& transaction_factory)
{
    NULL_USE(use_time_refinement);
    
    TBOX_ASSERT(dst_level);
    TBOX_ASSERT((next_coarser_level == -1) || hierarchy);
    
#ifdef HAMERS_DEBUG_CHECK_DIM_ASSERTIONS
    if (src_level)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *src_level);
    }
    if (hierarchy)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *hierarchy);
    }
#endif
    
    // Do we all agree on the destination box_level?
    if (src_level)
    {
        if (next_coarser_level >= 0)
        {
        }
    }
    
    d_schedule_created = true;
    
    boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory> trans_factory(
        transaction_factory);
    
    if (!trans_factory)
    {
        trans_factory.reset(new SAMRAI::xfer::StandardRefineTransactionFactory);
    }
    
    boost::shared_ptr<SAMRAI::xfer::PatchLevelFullFillPattern> fill_pattern(
        boost::make_shared<SAMRAI::xfer::PatchLevelFullFillPattern>());
    
    return boost::make_shared<ExtendedRefineSchedule>(
        fill_pattern,
        dst_level,
        src_level,
        next_coarser_level,
        hierarchy,
        d_refine_classes,
        trans_factory,
        patch_strategy,
        false);
}


/*
 **************************************************************************************************
 *
 * Create a communication schedule that copies data from the interiors of the old level and coarser
 * levels into the ghost cells and interior cells of the given new level.
 *
 **************************************************************************************************
 */
boost::shared_ptr<ExtendedRefineSchedule>
ExtendedRefineAlgorithm::createSchedule(
    const boost::shared_ptr<SAMRAI::xfer::PatchLevelFillPattern>& fill_pattern,
    const boost::shared_ptr<SAMRAI::hier::PatchLevel>& dst_level,
    const boost::shared_ptr<SAMRAI::hier::PatchLevel>& src_level,
    const int next_coarser_level,
    const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
    SAMRAI::xfer::RefinePatchStrategy* patch_strategy,
    bool use_time_refinement,
    const boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory>& transaction_factory)
{
    NULL_USE(use_time_refinement);
    
    TBOX_ASSERT(dst_level);
    TBOX_ASSERT((next_coarser_level == -1) || hierarchy);
    
#ifdef HAMERS_DEBUG_CHECK_DIM_ASSERTIONS
    if (src_level)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *src_level);
    }
    if (hierarchy)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY2(*dst_level, *hierarchy);
    }
#endif
    
    // Do we all agree on the destination box_level?
    if (src_level)
    {
       if (next_coarser_level >= 0)
       {
       }
    }
    
    d_schedule_created = true;
    
    boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory> trans_factory(
        transaction_factory);
    
    if (!trans_factory)
    {
        trans_factory.reset(new SAMRAI::xfer::StandardRefineTransactionFactory);
    }
    
    return boost::make_shared<ExtendedRefineSchedule>(
        fill_pattern,
        dst_level,
        src_level,
        next_coarser_level,
        hierarchy,
        d_refine_classes,
        trans_factory,
        patch_strategy,
        false);
}


/*
 **************************************************************************************************
 *
 * Reconfigure refine schedule to perform operations in this algorithm.
 *
 **************************************************************************************************
 */
bool
ExtendedRefineAlgorithm::checkConsistency(
    const boost::shared_ptr<ExtendedRefineSchedule>& schedule) const
{
    TBOX_ASSERT(schedule);
    return d_refine_classes->classesMatch(schedule->getEquivalenceClasses());
}


void ExtendedRefineAlgorithm::resetSchedule(
    const boost::shared_ptr<ExtendedRefineSchedule>& schedule) const
{
    TBOX_ASSERT(schedule);
    if (d_refine_classes->classesMatch(schedule->getEquivalenceClasses()))
    {
        schedule->reset(d_refine_classes);
    }
    else
    {
        TBOX_ERROR("ExtendedRefineAlgorithm::resetSchedule error..."
            << "\n Items in RefineClasses object passed to reset"
            << "\n routine does not match that in existing schedule."
            << std::endl);
    }
}


/*
 **************************************************************************************************
 *
 * Print refine algorithm data to the specified output stream.
 *
 **************************************************************************************************
 */
void
ExtendedRefineAlgorithm::printClassData(
   std::ostream& stream) const
{
    stream << "ExtendedRefineAlgorithm::printClassData()" << std::endl;
    stream << "----------------------------------------" << std::endl;
    d_refine_classes->printClassData(stream);
}
