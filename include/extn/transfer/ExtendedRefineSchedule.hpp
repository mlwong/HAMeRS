/*************************************************************************
 *
 * This file is modified from RefineSchedule.h of the SAMRAI version 3.11.2
 * distribution.  For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of the SAMRAI distribution.
 *
 ************************************************************************/

#ifndef EXTENDED_REFINE_SCHEDULE
#define EXTENDED_REFINE_SCHEDULE

#include "HAMeRS_config.hpp"

#include "SAMRAI/xfer/PatchLevelFillPattern.h"
#include "SAMRAI/xfer/RefineClasses.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/RefineTransactionFactory.h"
#include "SAMRAI/xfer/SingularityPatchStrategy.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Schedule.h"
#include "SAMRAI/tbox/Timer.h"

#include "boost/shared_ptr.hpp"
#include <iostream>

/*!
 * @brief Class ExtendedRefineSchedule performs the communication operations that refine data to,
 * copy data to, or fill physical boundary data on a destination patch level.
 *
 * Source data is copied into the provided scratch space for temporary processing.  The scratch space
 * must contain sufficient ghost cells to accommodate the stencil width of the given interpolation
 * operators and any physical boundary data that must be filled.  The scratch data is copied into
 * the destination space at the end of the process. The communication schedule is executed by calling
 * member function fillData().
 *
 * Each schedule object is typically created by a refine algorithm and represents the communication
 * dependencies for a particular configuration of the AMR hierarchy.  The communication schedule is
 * only valid for that particular configuration and must be regenerated when the AMR patch hierarchy
 * changes.  However, as long as the patch levels involved in the creation of the schedule remain
 * unchanged, the schedule may be used for multiple communication cycles.  For more information about
 * creating refine schedules, see the RefineAlgorithm header file.
 *
 * Some constructors accept the argument @c dst_level_fill_pattern.  This is a PatchLevelFillPattern
 * which controls which types of cells are filled and which are omitted from the filling process.
 * Concrete implementations of PatchLevelFillPattern are:
 * - @c PatchLevelFullFillPattern - Fill interior and ghost cells.
 * - @c PatchLevelInteriorFillPattern - Fill interior cells only.
 * - @c PatchLevelBorderFillPattern - Fill ghosts on level borders only.
 * - @c PatchLevelBorderAndInteriorFillPattern - Fill interior and
 *      ghosts on level borders.
 *
 * @see RefineAlgorithm
 * @see RefinePatchStrategy
 * @see RefineClasses
 */

class ExtendedRefineSchedule
{
    public:
        /*!
         * @brief Constructor that creates a refine schedule to copy data from the interiors of the
         * source patch data on the source level into the interiors and ghosts of destination patch
         * data on the destination level.
         *
         * The fill pattern supplied may restrict the data that will be copied.
         *
         * Only data on the intersection of the source and destination patch data will be copied; no
         * interpolation from coarser levels will be done.  The source and destination patch levels
         * must reside in the same index space.  However, the levels do not have to be in the same
         * AMR patch hierarchy.  Generally, this constructor is called by a RefineAlgorithm object.
         *
         * @param[in] dst_level_fill_pattern  Indicates which parts of the destination level to fill.
         * @param[in] dst_level               boost::shared_ptr to destination patch level.
         * @param[in] src_level               boost::shared_ptr to source patch level.
         * @param[in] refine_classes          boost::shared_ptr to structure containing patch data
         *                                    and operator information.  In general, this is
         *                                    constructed by the calling RefineAlgorithm object.
         * @param[in] transaction_factory     boost::shared_ptr to a factory object that will create
         *                                    data transactions.
         * @param[in] patch_strategy          Pointer to a refine patch strategy object that provides
         *                                    user-defined physical boundary filling operations.
         *                                    This pointer may be null, in which case no boundary
         *                                    filling operations will occur.  If your mesh has a
         *                                    singularity, the object this points to should have also
         *                                    inherited from SingularityPatchStrategy.
         * @param[in] use_time_interpolation  Boolean flag indicating whether to use time interpolation
         *                                    when setting data on the destination level. Default
         *                                    is no time interpolation.
         *
         * @pre dst_level
         * @pre src_level
         * @pre refine_classes
         * @pre transaction_factory
         * @pre dst_level->getDim() == src_level->getDim()
         * @pre dst_level->getGridGeometry()->getNumberOfBlockSingularities() == 0 || d_singularity_patch_strategy
         */
        ExtendedRefineSchedule(
            const boost::shared_ptr<SAMRAI::xfer::PatchLevelFillPattern>& dst_level_fill_pattern,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& dst_level,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& src_level,
            const boost::shared_ptr<SAMRAI::xfer::RefineClasses>& refine_classes,
            const boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory>& transaction_factory,
            SAMRAI::xfer::RefinePatchStrategy* patch_strategy,
            bool use_time_interpolation = false);
        
        /*!
         * @brief Constructor that creates a refine schedule to fill destination patch data on the
         * destination level from source patch data on the source level as well as interpolated data
         * from coarser levels.
         *
         * The fill pattern supplied may restrict the data that will be copied.
         *
         * Only data on the intersection of the source and destination patch data will be copied.
         * If portions of the destination level remain unfilled, then the algorithm recursively fills
         * those unfilled portions by interpolating source data from coarser levels in the AMR
         * hierarchy.  The source and destination patch levels must reside in the same index space.
         * However, the levels do not have to be in the same AMR patch hierarchy.  In general, this
         * constructor is called by a RefineAlgorithm object.
         *
         * @param[in] dst_level_fill_pattern  Indicates which parts of the destination level to fill.
         * @param[in] dst_level               boost::shared_ptr to destination patch level.  This
         *                                    level may be a level on the hierarchy or a coarsened
         *                                    version.
         * @param[in] src_level               boost::shared_ptr to source patch level; must be in
         *                                    same index space as destination level.  This pointer
         *                                    may be null, in which case the destination level will
         *                                    be filled only using data interpolated from coarser
         *                                    levels in the AMR hierarchy.
         * @param[in] next_coarser_level      Level number of next coarser level in AMR patch hierarchy
         *                                    relative to the destination level.  Note that when the
         *                                    destination level has number zero (i.e., the coarsest
         *                                    level), this value should be less than zero.
         * @param[in] hierarchy               boost::shared_ptr to patch hierarchy.  This pointer may
         *                                    be null only if the next_coarser_level value is < 0,
         *                                    indicating that there is no level in the hierarchy
         *                                    coarser than the destination level.
         * @param[in] refine_classes          boost::shared_ptr to structure containing patch data
         *                                    and operator information.  In general, this is
         *                                    constructed by the calling RefineAlgorithm object.
         * @param[in] transaction_factory     boost::shared_ptr to a factory object that will create
         *                                    data transactions.
         * @param[in] patch_strategy          Pointer to a refine patch strategy object that provides
         *                                    user-defined physical boundary filling operations.
         *                                    This pointer may be null, in which case no boundary
         *                                    filling or user-defined refine operations will occur.
         *                                    If your mesh has a singularity, the object this points
         *                                    to should have also inherited from
         *                                    SingularityPatchStrategy.
         * @param[in] use_time_refinement     Boolean flag indicating whether to use time interpolation
         *                                    when setting data on the destination level.  Default is
         *                                    no time interpolation.
         *
         * @pre dst_level
         * @pre (next_coarser_level == -1) || hierarchy
         * @pre refine_classes
         * @pre !src_level || (dst_level->getDim() == src_level.getDim())
         * @pre !hierarchy || (dst_level->getDim() == hierarchy.getDim())
         * @pre dst_level->getGridGeometry()->getNumberOfBlockSingularities() == 0 || d_singularity_patch_strategy
         */
        ExtendedRefineSchedule(
            const boost::shared_ptr<SAMRAI::xfer::PatchLevelFillPattern>& dst_level_fill_pattern,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& dst_level,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& src_level,
            int next_coarser_level,
            const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
            const boost::shared_ptr<SAMRAI::xfer::RefineClasses>& refine_classes,
            const boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory>& transaction_factory,
            SAMRAI::xfer::RefinePatchStrategy* patch_strategy,
            bool use_time_refinement = false);
        
        /*!
         * Destructor for the schedule releases all internal storage.
         */
        ~ExtendedRefineSchedule();
        
        /*!
         * @brief Reset this refine schedule to perform data transfers asssociated with refine class
         * items in function argument.
         *
         * In general, this function is called by a RefineAlgorithm object, which first checks that
         * the refine_classes parameter is in a state consistent with the ExtendedRefineSchedule
         * object.
         *
         * @param[in] refine_classes  boost::shared_ptr to structure containing patch data and operator
         *                            information.  In general, this is constructed by the calling
         *                            RefineAlgorithm object.  This pointer must be non-null.
         *
         * @pre refine_classes
         */
        void
        reset(
            const boost::shared_ptr<SAMRAI::xfer::RefineClasses>& refine_classes);
        
        /*!
         * @brief Execute the stored communication schedule and perform the data movement.
         *
         * @param[in] fill_time                  Time for filling operation.
         * @param[in] do_physical_boundary_fill  Boolean flag that can be used to bypass the physical
         *                                       boundary data filling operations on the destination
         *                                       level.  The default value is true indicating that
         *                                       boundary data will be filled (assuming a non-null
         *                                       refine patch strategy pointer was passed to the
         *                                       createSchedule() function.  Note that even when the
         *                                       value is false, boundary routines may be called on
         *                                       levels coarser than the destination level if such
         *                                       data is needed for proper interpolation.
         */
        void
        fillData(
            double fill_time,
            bool do_physical_boundary_fill = true) const;
        
        /*!
         * @brief Begin the stored communication schedule and perform the data movement. The
         * communication of data transfer in the same level is non-blocking. This method must be
         * followed by a call to <TT>finalizeNonBlockingfillData</TT> in order to complete the
         * communication.
         *
         * @param[in] fill_time                  Time for filling operation.
         * @param[in] do_physical_boundary_fill  Boolean flag that can be used to bypass the physical
         *                                       boundary data filling operations on the destination
         *                                       level.  The default value is true indicating that
         *                                       boundary data will be filled (assuming a non-null
         *                                       refine patch strategy pointer was passed to the
         *                                       createSchedule() function.  Note that even when the
         *                                       value is false, boundary routines may be called on
         *                                       levels coarser than the destination level if such
         *                                       data is needed for proper interpolation.
         */
        void
        beginNonBlockingfillData(
            double fill_time,
            bool do_physical_boundary_fill = true);
        
        /*!
         * @brief Finish the stored communication schedule and perform the data movement. This method
         * must follow <TT>finalizeNonBlockingfillData</TT> in order to complete the communication.
         */
        void
        finalizeNonBlockingfillData() const;
        
        /*!
         * @brief Return refine equivalence classes.
         *
         * The equivalence class information is used in schedule classes.
         */
        const boost::shared_ptr<SAMRAI::xfer::RefineClasses>&
        getEquivalenceClasses() const
        {
            return d_refine_classes;
        }
        
        /*!
         * @brief Set whether to unpack messages in a deterministic order.
         *
         * By default message unpacking is ordered by receive time, which is not deterministic.  If
         * your results are dependent on unpack ordering and you want deterministic results, set
         * this flag to true.
         *
         * @param [in] flag
         */
        void
        setDeterministicUnpackOrderingFlag(
            bool flag);
        
        /*!
         * @brief Print the refine schedule data to the specified data stream.
         *
         * @param[out] stream Output data stream.
         */
        void
        printClassData(
            std::ostream& stream) const;
        
    private:
        /*
         * Static integer constant describing the largest possible ghost cell width.
         */
        static const int BIG_GHOST_CELL_WIDTH = 10;
        
        ExtendedRefineSchedule(
            const ExtendedRefineSchedule&);                   // not implemented
        ExtendedRefineSchedule&
        operator = (
            const ExtendedRefineSchedule&);                   // not implemented
        
        /*!
         * @brief Allocate static timers.
         *
         * Only called by StartupShutdownManager.
         */
        static void
        initializeCallback();
        
        /*!
         * @brief Free static timers.
         *
         * Only called by StartupShutdownManager.
         */
        static void
        finalizeCallback();
        
        //! @brief Mapping from a (potentially remote) Box to a set of neighbors.
        typedef std::map<SAMRAI::hier::Box, SAMRAI::hier::BoxContainer, SAMRAI::hier::Box::id_less>
            FullNeighborhoodSet;
        
        /*!
         * @brief This private constructor creates a communication schedule that fills the destination
         * level interior as well as ghost regions equal to the maximum stencil width for refinement
         * operations.
         *
         * This constructor is used by the refine schedule algorithm during the recursive schedule
         * generation process.
         *
         * @param[out] errf Error flag
         * @param[in] dst_level               A temporary level that is used during interpolation.
         * @param[in] src_level               A level from the hierarchy that is of the same resolution
         *                                    as dst_level
         * @param[in] next_coarser_level      Level number of next coarser level in AMR patch hierarchy
         *                                    relative to the destination level.  Note that when the
         *                                    destination level has number zero (i.e., the coarsest
         *                                    level), this value should be less than zero.
         * @param[in] hierarchy               boost::shared_ptr to patch hierarchy.
         * @param[in] dst_to_src
         * @param[in] src_growth_to_nest_dst  The minimum amount that src_level has to grow in order
         *                                    to nest dst.
         * @param[in] refine_classes          Holds refine equivalence classes to be used by this
         *                                    schedule.
         * @param[in] transaction_factory     boost::shared_ptr to a factory object that will create
         *                                    data transactions.
         * @param[in] patch_strategy          boost::shared_ptr to a refine patch strategy object
         *                                    that provides user-defined physical boundary filling
         *                                    operations.  This pointer may be null, in which case
         *                                    no boundary filling or user-defined refine operations
         *                                    will occur.
         * @param[in] top_refine_schedule
         *
         * @pre dst_level
         * @pre src_level
         * @pre (next_coarser_level == -1) || hierarchy
         * @pre dst_to_src.hasTranspose()
         * @pre refine_classes
         * @pre dst_level->getDim() == src_level.getDim()
         * @pre !hierarchy || (dst_level->getDim() == hierarchy.getDim())
         * @pre dst_to_src.getBase() == *dst_level->getBoxLevel()
         * @pre src_to_dst.getHead() == *dst_level->getBoxLevel()
         */
        ExtendedRefineSchedule(
            int& errf,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& dst_level,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& src_level,
            int next_coarser_level,
            const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
            const SAMRAI::hier::Connector& dst_to_src,
            const SAMRAI::hier::IntVector& src_growth_to_nest_dst,
            const boost::shared_ptr<SAMRAI::xfer::RefineClasses>& refine_classes,
            const boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory>& transaction_factory,
            SAMRAI::xfer::RefinePatchStrategy* patch_strategy,
            const ExtendedRefineSchedule* top_refine_schedule);
        
        /*!
         * @brief Read static data from input database.
         */
        void
        getFromInput();
        
        /*!
         * @brief Finish the schedule construction for the two constructors that take a hierarchy
         * as an argument.
         *
         * The hierarchy gives the possibility of recursion to get data from coarser levels.
         *
         * Returns an error flag:
         * - 0: No error.
         * - 1: Internal error due to unfilled boxes without coarser level.
         * - 2: Internal error due to unfilled boxes without hierarchy.
         *
         * @param[in] next_coarser_ln         Level number of the level coarser than the destination
         *                                    level
         * @param[in] hierarchy               A patch hierarchy to be used to provide coarser levels
         *                                    that may be used for interpolation.
         * @param[in] src_growth_to_nest_dst  The minimum amount that the source level has to grow
         *                                    in order to nest the destination level.
         * @param[in] dst_to_fill             Connector describing the boxes to fill for each
         *                                    destination patch.
         * @param[in] src_owner_dst_to_fill   A BoxNeighborhoodCollection that maps each local box
         *                                    on the source level to a collection of boxes that
         *                                    indicates what parts of the destination fill boxes can
         *                                    be filled by that source box.
         * @param[in] use_time_interpolation  Boolean flag indicating whether to use time interpolation
         *                                    when setting data on the destination level.
         * @param[in] skip_generate_schedule  If true, then the generation of transactions to
         *                                    communicate from source level to destination level will
         *                                    be skipped.
         *
         * @pre d_dst_to_src
         * @pre d_dst_to_src->hasTranspose()
         * @pre (next_coarser_ln == -1) || hierarchy
         * @pre !d_src_level || d_dst_to_src->isFinalized()
         * @return error flag
         */
        int
        finishScheduleConstruction(
            int next_coarser_ln,
            const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
            const SAMRAI::hier::IntVector& src_growth_to_nest_dst,
            const SAMRAI::hier::Connector& dst_to_fill,
            const SAMRAI::hier::BoxNeighborhoodCollection& src_owner_dst_to_fill,
            bool use_time_interpolation,
            bool skip_generate_schedule = false);
        
        /*!
         * @brief Allocate scratch space on the specified level and store the allocated patch data
         * indices in the component selector for later deallocation.
         *
         * @param[out] allocate_vector  Component selector that will store the allocated patch data
         *                              indices.
         * @param[in,out] level
         * @param[in] fill_time         Simulation time for filling operation.
         *
         * @pre level
         */
        void
        allocateScratchSpace(
            SAMRAI::hier::ComponentSelector& allocate_vector,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& level,
            double fill_time) const;
        
        void
        allocateDestinationSpace(
            SAMRAI::hier::ComponentSelector& allocate_vector,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& level,
            double fill_time) const;
        
        /*!
         * @brief Recursively fill the destination level with data at the given time.
         *
         * @param[in]  fill_time                  Simulation time when the fill takes place
         * @param[in]  do_physical_boundary_fill  Indicates whether to call user-supplied boundary
         *                                        filling routines regardless of whether this is
         *                                        needed based on ghost cell width of destination
         *                                        data or stencil width of some interpolation
         *                                        operator.
         */
        void
        recursiveFill(
            double fill_time,
            bool do_physical_boundary_fill) const;
        
        /*!
         * @brief Begin the recursive filling process with non-blocking communication for data
         * transfer in the same level of refinement.
         *
         * This method must be followed by a call to <TT>finalizeNonBlockingRecursiveFill</TT> in
         * order to complete the filling process for data on the same level.
         * 
         * @param[in]  fill_time                  Simulation time when the fill takes place
         * @param[in]  do_physical_boundary_fill  Indicates whether to call user-supplied boundary
         *                                        filling routines regardless of whether this is
         *                                        needed based on ghost cell width of destination
         *                                        data or stencil width of some interpolation
         *                                        operator.
         *
         */
        void
        beginNonBlockingRecursiveFill(
            double fill_time,
            bool do_physical_boundary_fill);
        
        /*!
         * @brief Finish the recursive filling process with non-blocking communication for data
         * transfer in the same level of refinement.
         *
         * This method must follow <TT>beginNonBlockingRecursiveFill</TT> and filling process is
         * completed when this method returns.
         * 
         */
        void
        finalizeNonBlockingRecursiveFill() const;
        
        /*!
         * @brief Fill the physical boundaries for each patch on d_dst_level.
         *
         * @param[in] fill_time  Simulation time when the fill takes place
         *
         * @pre d_dst_level
         */
        void
        fillPhysicalBoundaries(
            double fill_time) const;
     
        void
        fillSingularityBoundaries(
            double fill_time) const;
        
        /*!
         * @brief Copy the scratch space into the destination space in d_dst_level.
         *
         * If the scratch and destination patch data components are the same, then no copying is
         * performed.
         *
         * @pre d_dst_level
         */
        void
        copyScratchToDestination() const;
        
        /*!
         * @brief Refine scratch data between coarse and fine patch levels.
         *
         * @param[in] fine_level          Fine level to receive interpolated data
         * @param[in] coarse_level        Coarse level source of interpolation
         * @param[in] coarse_to_fine      Connector coarse to fine
         * @param[in] coarse_to_unfilled  Connector coarse to level representing boxes that need to
         *                                be filled.
         * @param[in] overlaps
         */
        void
        refineScratchData(
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& fine_level,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& coarse_level,
            const SAMRAI::hier::Connector& coarse_to_fine,
            const SAMRAI::hier::Connector& coarse_to_unfilled,
            const std::vector<std::vector<boost::shared_ptr<SAMRAI::hier::BoxOverlap> > >&
            overlaps) const;
        
        /*!
         * @brief Compute and store the BoxOverlaps that will be needed by refineScratchData().
         */
        void
        computeRefineOverlaps(
            std::vector<std::vector<boost::shared_ptr<SAMRAI::hier::BoxOverlap> > >& overlaps,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& fine_level,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& coarse_level,
            const SAMRAI::hier::Connector& coarse_to_fine,
            const SAMRAI::hier::Connector& coarse_to_unfilled);
        
        /*!
         * @brief Constructs the transactions for all communication and copying of patch data.
         *
         * The resulting transactions will only fill the regions of intersection
         * between the fill level, which is "head" level that the Connector
         * dst_to_fill points to, and the source level.  The remaining
         * box regions are added to unfilled_box_level.
         *
         * @param[out] unfilled_box_level        The parts of the fill level that cannot be filled
         *                                       from the source level are added here.
         * @param[out] dst_to_unfilled           Connector from dst_level to unfilled_box_level.
         * @param[out] unfilled_encon_box_level  The parts of the fill level at enhanced connectivity
         *                                       block boundaries that cannot be filled from the
         *                                       source level.
         * @param[out] encon_to_unfilled_encon   Connector from level representing enhanced connectivity
         *                                       on the destination level to unfilled_encon_box_level
         * @param[in] dst_to_fill                Connector between dst_level and a level representing
         *                                       the boxes that need to be filled.
         * @param[in] src_owner_dst_to_fill      A BoxNeighborhoodCollection that maps each local
         *                                       box on the source level to a collection of boxes
         *                                       that indicates what parts of the fill can be filled
         *                                       by that source box.
         * @param[in] use_time_interpolation     Boolean flag indicating whether to use time
         *                                       interpolation when setting data on the destination
         *                                       level.
         *
         * @param[in] create_transactions
         *
         * @pre d_dst_to_src
         * @pre d_dst_to_src->hasTranspose()
         */
        void
        generateCommunicationSchedule(
            boost::shared_ptr<SAMRAI::hier::BoxLevel>& unfilled_box_level,
            boost::shared_ptr<SAMRAI::hier::Connector>& dst_to_unfilled,
            boost::shared_ptr<SAMRAI::hier::BoxLevel>& unfilled_encon_box_level,
            boost::shared_ptr<SAMRAI::hier::Connector>& encon_to_unfilled_encon,
            const SAMRAI::hier::Connector& dst_to_fill,
            const SAMRAI::hier::BoxNeighborhoodCollection& src_owner_dst_to_fill,
            const bool use_time_interpolation,
            const bool create_transactions);
        
        /*!
         * @brief Compute boxes that need to be filled and data associated with them.
         *
         * fill_box_level will be filled with boxes representing all of the regions intended to be
         * filled by the schedule.  It will include the boxes of dst_box_level grown by ill_gcw, but
         * then can be restricted based on the PatchLevelFillPattern given to the schedule constructor.
         * This method sets up this fill level, as well as a connector from the destination level to
         * the fill level, and also provides a mapping from each source box to the boxes they can
         * fill on the fill level.
         *
         * @param[out] fill_box_level         Will contain all boxes that need to be filled.
         * @param[out] dst_to_fill            Connector from dst_level to fill_box_level.
         * @param[out] src_owner_dst_to_fill  A BoxNeighborhoodCollection that maps each local box
         *                                    on the source level to a collection of boxes that
         *                                    indicates what parts of fill_box_level can be filled
         *                                    by that source box.
         *
         * @pre (d_dst_level->getDim() == dst_box_level.getDim()) &&
         *      (d_dst_level->getDim() == fill_ghost_width.getDim())
         *
         * @pre d_dst_to_src
         * @pre d_dst_to_src->hasTranspose()
         */
        void
        setDefaultFillBoxLevel(
            boost::shared_ptr<SAMRAI::hier::BoxLevel>& fill_box_level,
            boost::shared_ptr<SAMRAI::hier::Connector>& dst_to_fill,
            SAMRAI::hier::BoxNeighborhoodCollection& src_owner_dst_to_fill);
        
        /*
         * @brief Set up level to represent ghost regions at enhanced connectivity block boundaries.
         *
         * @param[in] fill_gcw  Width to extend across the block boundary
         */
        void
        createEnconLevel(
            const SAMRAI::hier::IntVector& fill_gcw);
        
        /*
         * @brief Find the fill boxes that are at enhanced connectivity.
         *
         * Given a list representing fill boxes, determine the portion of those boxes that lie across
         * any enhanced connectivity boundary from the block specified by dst_block_id.
         *
         * @param[out]  encon_fill_boxes
         * @param[in]   fill_boxes_list
         * @param[in]   dst_block_id
         *
         * @pre encon_fill_boxes.empty()
         */
        void
        findEnconFillBoxes(
            SAMRAI::hier::BoxContainer& encon_fill_boxes,
            const SAMRAI::hier::BoxContainer& fill_boxes_list,
            const SAMRAI::hier::BlockId& dst_block_id);
        
        /*
         * @brief Find unfilled boxes at enhanced connectivity.
         *
         * Determine which portion of encon_fill_boxes cannot be filled from a source level.  Those
         * unfilled boxes are added to level_encon_unfilled_boxes, and edges are added to
         * encon_to_unfilled_encon_nbrhood_set.
         *
         * The source level is the head level from the connector d_dst_to_src.
         *
         * @param[out]  unfilled_encon_box_level    encon unfilled boxes for the dst level
         * @param[out]  encon_to_unfilled_encon     connector from d_encon_level to
         *                                          level_encon_unfilled_boxes
         * @param[in,out]  last_unfilled_local_id   a unique LocalId not already used in
         *                                          level_encon_unfilled_boxes
         * @param[in]  dst_box                      The destination box
         * @param[in]  encon_fill_boxes
         *
         * @pre d_dst_to_src
         */
        void
        findEnconUnfilledBoxes(
            const boost::shared_ptr<SAMRAI::hier::BoxLevel>& unfilled_encon_box_level,
            const boost::shared_ptr<SAMRAI::hier::Connector>& encon_to_unfilled_encon,
            SAMRAI::hier::LocalId& last_unfilled_local_id,
            const SAMRAI::hier::Box& dst_box,
            const SAMRAI::hier::BoxContainer& encon_fill_boxes);
        
        /*
         * @brief Create schedule for filling unfilled boxes at enhanced connectivity.
         *
         * @param[in]  hierarchy                The patch hierarchy
         * @param[in]  hiercoarse_level         Level on hierarchy one level coarser than the
         *                                      destination level
         * @param[in]  src_growth_to_nest_dst
         * @param[in]  encon_to_unfilled_encon
         */
        void
        createEnconFillSchedule(
            const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& hiercoarse_level,
            const SAMRAI::hier::IntVector& src_growth_to_nest_dst,
            const SAMRAI::hier::Connector& encon_to_unfilled_encon);
        
        /*!
         * @brief Communicate dst_to_fill info to the src owners when the owners would otherwise be
         * unable to compute the info.
         *
         * @param[out] src_owner_dst_to_fill  A BoxNeighborhoodCollection that maps each local box
         *                                    on the source level to a collection of boxes that
         *                                    indicates what parts of fill_box_level can be filled
         *                                    by that source box.
         * @param[in] dst_to_fill             Mapping from the dst_level to boxes it needs need to
         *                                    have filled.
         *
         * @pre d_dst_to_src
         * @pre d_dst_to_src->hasTranspose()
         */
        void
        communicateFillBoxes(
            SAMRAI::hier::BoxNeighborhoodCollection& src_owner_dst_to_fill,
            const SAMRAI::hier::Connector& dst_to_fill);
        
        /*!
         * @brief Shear off parts of unfilled boxes that lie outside non-periodic domain boundaries
         * and update an overlap Connector based on the change.
         *
         * @param[in,out] unfilled
         *
         * @param[in,out] dst_to_unfilled
         *
         * @param[in] hierarchy
         */
        void
        shearUnfilledBoxesOutsideNonperiodicBoundaries(
            SAMRAI::hier::BoxLevel& unfilled,
            SAMRAI::hier::Connector& dst_to_unfilled,
            const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy);
        
        /*!
         * @brief Make an unfilled box level consisting of node-centered boxes.
         *
         * When this is called, there already exists a cell-centered unfilled box level, which is
         * the first argument passed in here. This method creates a node-centered unfilled box level,
         * consisting of the boxes in the cell-centered unfilled level converted to a node centering,
         * with any nodes that can be filled by this schedule's source level removed.
         *
         * @param[in] unfilled_box_level  Cell-centered unfilled level
         * @param[in] dst_to_unfiled      Connector from destination to unfilled
         */
        void
        makeNodeCenteredUnfilledBoxLevel(
            const SAMRAI::hier::BoxLevel& unfilled_box_level,
            const SAMRAI::hier::Connector& dst_to_unfilled);
        
        /*!
         * @brief Set up a coarse interpolation BoxLevel and related data.
         *
         * Also sets up dst_to_coarse_interp, coarse_interp_to_dst and coarse_interp_to_unfilled.
         *
         * @param[out] coarse_interp_box_level
         *
         * @param[out] dst_to_coarse_interp
         *
         * @param[out] coarse_interp_to_unfilled
         *
         * @param[in] hiercoarse_box_level The BoxLevel on the hierarchy at the resolution that
         * coarse_interp_box_level is to have.
         *
         * @param[in] dst_to_unfilled
         */
        void
        setupCoarseInterpBoxLevel(
            boost::shared_ptr<SAMRAI::hier::BoxLevel>& coarse_interp_box_level,
            boost::shared_ptr<SAMRAI::hier::Connector>& dst_to_coarse_interp,
            boost::shared_ptr<SAMRAI::hier::Connector>& coarse_interp_to_unfilled,
            const SAMRAI::hier::BoxLevel& hiercoarse_box_level,
            const SAMRAI::hier::Connector& dst_to_unfilled);
        
        /*!
         * @brief Create a coarse interpolation PatchLevel and compute the Connectors between the
         * coarse interpolation level and the hiercoarse level.
         *
         * @param[out] coarse_interp_level
         *
         * @param[in,out] coarse_interp_box_level This method will add periodic images to
         * coarse_interp_box_level, if needed.
         *
         * @param[in] coarse_interp_to_hiercoarse
         *
         * @param[in] next_coarser_ln Level number of hiercoarse (the coarser level on the hierarchy)
         *
         * @param[in] hierarchy
         *
         * @param[in] dst_to_src
         *
         * @param[in] dst_to_coarse_interp
         *
         * @param[in] dst_level
         *
         * @pre dst_to_src.hasTranspose()
         * @pre dst_to_coarse_interp.hasTranspose()
         */
        void
        createCoarseInterpPatchLevel(
            boost::shared_ptr<SAMRAI::hier::PatchLevel>& coarse_interp_level,
            boost::shared_ptr<SAMRAI::hier::BoxLevel>& coarse_interp_box_level,
            boost::shared_ptr<SAMRAI::hier::Connector>& coarse_interp_to_hiercoarse,
            const int next_coarser_ln,
            const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
            const SAMRAI::hier::Connector& dst_to_src,
            const SAMRAI::hier::Connector& dst_to_coarse_interp,
            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& dst_level);
        
        /*!
         * @brief Check that the Connectors between the coarse interpolation and hiercoarse levels
         * are transposes and that the coarse interpolation BoxLevel sufficiently nests inside the
         * hiercoarse.
         */
        void
        sanityCheckCoarseInterpAndHiercoarseLevels(
            const int next_coarser_ln,
            const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
            const SAMRAI::hier::Connector& coarse_interp_to_hiercoarse);
        
        /*!
         * @brief Get whether there is data living on patch borders.
         */
        bool
        getDataOnPatchBorderFlag() const;
        
        /*!
         * @brief Get the mininum width required for the overlap Connector between destination and
         * source.
         */
        SAMRAI::hier::IntVector
        getMinConnectorWidth() const;
        
        /*!
         * @brief Get the maximum ghost cell width of all destination patch data components.
         */
        SAMRAI::hier::IntVector
        getMaxDestinationGhosts() const;
        
        /*!
         * @brief Get the maximum ghost cell width of all scratch patch data components.
         */
        SAMRAI::hier::IntVector
        getMaxScratchGhosts() const;
        
        /*!
         * @brief Get the maximum ghost cell width required for all stencils.
         */
        SAMRAI::hier::IntVector
        getMaxStencilGhosts() const;
        
        /*!
         * @brief Function that constructs schedule transactions between one source box and one
         * destination box.
         *
         * Transactions will move data on the intersection of the source and destination boxes with
         * the fill boxes.
         *
         * @param[in] num_nbrs
         * @param[in] nbrs_begin
         * @param[in] nbrs_end
         * @param[in] dst_box  Box from a destination patch.
         * @param[in] src_box  Box from a source patch.
         * @param[in] use_time_interpolation
         *
         * @pre d_dst_level
         * @pre d_src_level
         * @pre !dst_box.isPeriodicImage()
         */
        void
        constructScheduleTransactions(
            int num_nbrs,
            SAMRAI::hier::BoxNeighborhoodCollection::ConstNeighborIterator& nbrs_begin,
            SAMRAI::hier::BoxNeighborhoodCollection::ConstNeighborIterator& nbrs_end,
            const SAMRAI::hier::Box& dst_box,
            const SAMRAI::hier::Box& src_box,
            const bool use_time_interpolation);
        
        /*!
         * @brief Reorder the neighborhood sets from a src_to_dst Connector so they can be used in
         * schedule generation.
         *
         * First, this puts the neighborhood set data in src_to_dst into dst-major order so the src
         * owners can easily loop through the dst-src edges in the same order that dst owners see
         * them.  Transactions must have the same order on the sending and receiving processors.
         *
         * Section, it shifts periodic image dst boxes back to the zero-shift position, and applies
         * a similar shift to src boxes so that the overlap is unchanged.  The
         * constructScheduleTransactions method requires all shifts to be absorbed in the src box.
         *
         * The reordered neighboorhood sets are added to the output parameter.
         *
         * @param[out] full_inverted_edges
         *
         * @pre d_dst_to_src
         * @pre d_dst_to_src->hasTranspose()
         */
        void
        reorderNeighborhoodSetsByDstNodes(
            FullNeighborhoodSet& full_inverted_edges) const;
        
        /*!
         * @brief Cache local copies of hierarchy information and compute necessary ghost and stencil
         * widths.
         *
         * This is called by every ExtendedRefineSchedule constructor.
         */
        void
        initializeDomainAndGhostInformation();
        
        /*!
         * @brief Utility function to set up local copies of refine items.
         *
         * An array of refine data items obtained from the CoarsenClasses object is stored locally
         * here to facilitate interaction with transactions.
         *
         * @param[in] refine_classes
         */
        void
        setRefineItems(
            const boost::shared_ptr<SAMRAI::xfer::RefineClasses>& refine_classes);
        
        /*
         * @brief Utility function to clear local copies of refine items.
         */
        void
        clearRefineItems();
        
        /*!
         * @brief Utility function to check refine items to see whether their patch data components
         * have sufficient ghost width to handle user-defined interpolation operations.
         *
         * Specifically scratch data ghost cell widths must be at least as large as the stencil of
         * those user-defined operations.
         *
         * If any of the tested ghost cell widths are insufficient, an error will occur with a
         * descriptive message.
         */
        void
        initialCheckRefineClassItems() const;
        
        /*!
         * Structures that store refine data items.
         */
        boost::shared_ptr<SAMRAI::xfer::RefineClasses> d_refine_classes;
        
        /*!
         * @brief number of refine data items
         */
        size_t d_number_refine_items;
        
        /*!
         * @brief used as array to store copy of refine data items.
         */
        const SAMRAI::xfer::RefineClasses::Data** d_refine_items;
        
        /*!
         * @brief boost::shared_ptr to the destination patch level.
         */
        boost::shared_ptr<SAMRAI::hier::PatchLevel> d_dst_level;
        
        /*!
         * @brief boost::shared_ptr to the source patch level.
         */
        boost::shared_ptr<SAMRAI::hier::PatchLevel> d_src_level;
        
        /*!
         * @brief Object supporting interface to user-defined boundary filling and spatial data
         * interpolation operations.
         */
        SAMRAI::xfer::RefinePatchStrategy* d_refine_patch_strategy;
        
        /*!
         * @brief Object supporting interface to user-defined ghost data filling at block
         * singularities.
         */
        SAMRAI::xfer::SingularityPatchStrategy* d_singularity_patch_strategy;
        
        /*!
         * @brief Factory object used to create data transactions when schedule is constructed.
         */
        boost::shared_ptr<SAMRAI::xfer::RefineTransactionFactory> d_transaction_factory;
        
        /*!
         * @brief  Whether there is data on patch borders.
         */
        bool d_data_on_patch_border_flag;
        
        /*!
         * @brief  maximum stencil width.
         */
        SAMRAI::hier::IntVector d_max_stencil_width;
        
        /*!
         * @brief maximum scratch ghost cell widths.
         */
        SAMRAI::hier::IntVector d_max_scratch_gcw;
        
        /*!
         * @brief Width of ghost cell region to fill passed to user supplied physical boundary
         * condition routine.
         */
        SAMRAI::hier::IntVector d_boundary_fill_ghost_width;
        
        /*!
         * @brief Flag indicating whether user's physical boundary data filling routine should be
         * forced at last step of level filling process.
         *
         * This flag is true when doing recursive filling, because the ghost data may be needed by
         * finer levels (regardless of whether the user requested ghost boundary filling).  This
         * variable is set in the constructors, which knows whether the object is being constructed
         * for recursive filling.
         *
         * For efficiency, we only force boundary filling when, during object construction, we
         * determine that the ghost cells do exist.
         */
        bool d_force_boundary_fill;
        
        /*!
         * @brief Boolean flag indicating whether physical domain can be represented as a single
         * box region.
         */
        std::vector<bool> d_domain_is_one_box;
        
        /*!
         * @brief Number of non-zero entries in periodic shift vector.
         */
        int d_num_periodic_directions;
        
        /*!
         * @brief the periodic shift vector.
         */
        SAMRAI::hier::IntVector d_periodic_shift;
        
        /*!
         * @brief  Level-to-level communication schedule between the source and destination.
         *
         * d_coarse_priority_level_schedule handles the situation where coarse data should take
         * precedence at coarse-fine boundaries for data types holding values at patch boundaries
         * but which are considered interior values.  d_fine_priority_level_schedule handles the
         * situation where fine data should take precedence.
         */
        boost::shared_ptr<SAMRAI::tbox::Schedule> d_coarse_priority_level_schedule;
        
        /*!
         * @brief Level-to-level communication schedule between the source and destination.
         *
         * d_coarse_priority_level_schedule handles the situation where coarse data should take
         * precedence at coarse-fine boundaries for data types holding values at patch boundaries
         * but which are considered interior values.  d_fine_priority_level_schedule handles the
         * situation where fine data should take precedence.
         */
        boost::shared_ptr<SAMRAI::tbox::Schedule> d_fine_priority_level_schedule;
        
        /*!
         * @brief The coarse interpolation level is an internal level created to hold data required
         * for interpolating into the fill boxes of the destination that could not be filled directly
         * from the source level.
         *
         * Once d_coarse_interp_level is filled (by executing d_coarse_interp_schedule) interpolating
         * data into the corresponding fill boxes of the destination is a local operation.
         *
         * This coarser level is filled by the d_coarse_interp_schedule.  If no coarser level data
         * is needed, then this pointer will be NULL.  Note that the coarse interpolation level may
         * not have the same mapping as the destination level.
         */
        boost::shared_ptr<SAMRAI::hier::PatchLevel> d_coarse_interp_level;
        
        /*!
         * @brief The coarse interpolation encon level is an internal level created to hold data
         * used for interpolating into unfilled boxes at enhanced connectivity block boundaries.
         *
         * d_coarse_interp_encon_level will be filled by d_coarse_interp_encon_schedule.  Once it
         * is filled, the interpolation of data to patches in d_encon_level will be a local
         * operation.
         */
        boost::shared_ptr<SAMRAI::hier::PatchLevel> d_coarse_interp_encon_level;
        
        /*!
         * @brief Schedule to recursively fill the coarse interpolation level using the next coarser
         * hierarchy level.
         *
         * This schedule describes how to fill the coarse interpolation level so that the coarse
         * data can be interpolated into the fine fill boxes on the destination.
         */
        boost::shared_ptr<ExtendedRefineSchedule> d_coarse_interp_schedule;
        
        /*!
         * @brief Schedule to recursively fill d_coarse_interp_encon_level using the next coarser
         * hierarchy level.
         *
         * This schedule fills d_coarse_interp_encon_level so that it can be used to interpolate
         * data onto d_encon_level in fill boxes that could not be filled from the source level.
         */
        boost::shared_ptr<ExtendedRefineSchedule> d_coarse_interp_encon_schedule;
        
        /*!
         * @brief Internal level representing ghost regions of destination patches at enhanced
         * connectivity block boundaries.
         *
         * When a destination patch touches an enhanced connectivity block boundary, a patch will
         * be created in the coordinate systems of its singularity neighbor blocks representing the
         * portion of the destination patch's ghost region that lies in those neighboring blocks.
         */
        boost::shared_ptr<SAMRAI::hier::PatchLevel> d_encon_level;
        
        /*!
         * @brief Intermediate destination level for interpolating ghost data at block boundaries.
         */
        boost::shared_ptr<SAMRAI::hier::PatchLevel> d_nbr_blk_fill_level;
        
        /*!
         * @brief Describes remaining unfilled boxes after attempting to fill from the source level.
         * These remaining boxes must be filled using a coarse interpolation schedule,
         * d_coarse_interp_schedule.
         */
        boost::shared_ptr<SAMRAI::hier::BoxLevel> d_unfilled_box_level;
        boost::shared_ptr<SAMRAI::hier::BoxLevel> d_unfilled_node_box_level;
        
        /*!
         * @brief Describes remaining unfilled boxes of d_encon_level after attempting to fill from
         * the source level.  These remaining boxes must be filled using a coarse interpolation
         * schedule, d_coarse_interp_encon_schedule.
         */
        boost::shared_ptr<SAMRAI::hier::BoxLevel> d_unfilled_encon_box_level;
        
        /*!
         * @brief Stores the BoxOverlaps needed by refineScratchData().
         */
        std::vector<std::vector<boost::shared_ptr<SAMRAI::hier::BoxOverlap> > >
        d_refine_overlaps;
        
        /*!
         * @brief Stores the BoxOverlaps needed by refineScratchData() for unfilled boxes at enhanced
         * connectivity.
         */
        std::vector<std::vector<boost::shared_ptr<SAMRAI::hier::BoxOverlap> > >
        d_encon_refine_overlaps;
        
        /*!
         * @brief Stores the overlaps needed to copy from d_nbr_blk_fill_level to the destination
         * at block boundaries.
         */
        std::vector<std::vector<boost::shared_ptr<SAMRAI::hier::BoxOverlap> > >
        d_nbr_blk_copy_overlaps;
        
        /*!
         * @brief Connector from the destination level to the coarse interpolation.
         */
        boost::shared_ptr<SAMRAI::hier::Connector> d_dst_to_coarse_interp;
        
        /*!
         * @brief Connector from d_encon_level to d_coarse_interp_encon_level.
         */
        boost::shared_ptr<SAMRAI::hier::Connector> d_encon_to_coarse_interp_encon;
        
        boost::shared_ptr<SAMRAI::hier::Connector> d_unfilled_to_unfilled_node;
        
        /*!
         * @brief Connector d_coarse_interp_level to d_unfilled_box_level.
         *
         * Cached for use during schedule filling.
         */
        boost::shared_ptr<SAMRAI::hier::Connector> d_coarse_interp_to_unfilled;
        
        boost::shared_ptr<SAMRAI::hier::Connector> d_coarse_interp_encon_to_unfilled_encon;
        
        /*!
         * @brief Connector from coarse interp level to d_nbr_blk_fill_level, used to connect coarse
         * and fine patches for interpolation at block boundaries.
         */
        boost::shared_ptr<SAMRAI::hier::Connector> d_coarse_interp_to_nbr_fill;
        
        boost::shared_ptr<SAMRAI::hier::Connector> d_dst_to_encon;
        boost::shared_ptr<SAMRAI::hier::Connector> d_encon_to_src;
        const SAMRAI::hier::Connector* d_dst_to_src;
        
        std::map<SAMRAI::hier::BoxId, SAMRAI::hier::IntVector> d_nbr_refine_ratio;
        std::map<SAMRAI::hier::BoxId, SAMRAI::hier::IntVector> d_encon_nbr_refine_ratio;
        
        //@{
        
        /*!
         * @name Data members used in constructScheduleTransactions().
         */
        
        /*!
         * @brief Array to store overlaps used for construction of transactions.
         *
         * This is declared in the class to make memory management more efficient since
         * constructScheduleTransactions is called many times.
         *
         * The size of the array is controlled by d_max_fill_boxes.
         */
        std::vector<boost::shared_ptr<SAMRAI::hier::BoxOverlap> > d_overlaps;
        
        /*!
         * @brief Source mask boxes used in construction of transactions.
         *
         * Like d_overlaps, this is declared in the class to make memory management more efficient.
         */
        SAMRAI::hier::BoxContainer d_src_masks;
        
        /*!
         * @brief The maximum number of fill boxes across all patches in the destination level.
         */
        int d_max_fill_boxes;
        
        //@}
        
        /*!
         * @brief PatchLevelFillPattern controlling what parts of the destination level can be filled.
         */
        boost::shared_ptr<SAMRAI::xfer::PatchLevelFillPattern> d_dst_level_fill_pattern;
        
        /*!
         * @brief Required fine Connector widths used in refining data from coarser levels.
         *
         * Maintained by the top level ExtendedRefineSchedule for use by the recursive schedules.
         * Unused when there is no hierarchy for recursion.
         */
        std::vector<SAMRAI::hier::IntVector> d_fine_connector_widths;
        
        /*!
         * @brief The top ExtendedRefineSchedule that led recursively to this one.
         *
         * The top ExtendedRefineSchedule contains some parameters to be shared by the recursive
         * ExtendedRefineSchedules.
         */
        const ExtendedRefineSchedule* d_top_refine_schedule;
        
        /*!
         * @brief Flag indicating whether user's physical boundary data filling routine should be
         * forced during the non-blocking data filling process.
         */
        bool d_do_physical_boundary_fill;
        
        /*!
         * @brief Time for non-blocking filling operation.
         */
        double d_fill_time;
        
        /*!
         * @brief Component selctors for non-blocking filling operation.
         */
        SAMRAI::hier::ComponentSelector d_allocate_vector;
        SAMRAI::hier::ComponentSelector d_encon_allocate_vector;
        SAMRAI::hier::ComponentSelector d_nbr_fill_scratch_vector;
        SAMRAI::hier::ComponentSelector d_nbr_fill_dst_vector;
        
        /*!
         * @brief Shared debug checking flag.
         */
        static bool s_extra_debug;
        
        /*!
         * @brief Flag that turns on barrier calls for use in performance analysis.
         */
        static bool s_barrier_and_time;
        
        /*!
         * @brief Flag indicating if any ExtendedRefineSchedule has read the input database for
         * static data.
         */
        static bool s_read_static_input;
        
        //@{
        /*!
         * @name Timer objects for performance measurement.
         */
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_refine_schedule;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_fill_data;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_fill_data_nonrecursive;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_fill_data_recursive;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_fill_physical_boundaries;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_fill_singularity_boundaries;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_refine_scratch_data;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_finish_sched_const;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_finish_sched_const_recurse;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_gen_comm_sched;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_shear;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_get_global_box_count;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_coarse_shear;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_setup_coarse_interp_box_level;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_bridge_coarse_interp_hiercoarse;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_bridge_dst_hiercoarse;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_invert_edges;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_construct_send_trans;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_construct_recv_trans;
        
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_coarse_priority_level_schedule_communicate;
        static boost::shared_ptr<SAMRAI::tbox::Timer> t_fine_priority_level_schedule_communicate;
        
        //@}
        
        static SAMRAI::tbox::StartupShutdownManager::Handler
            s_initialize_finalize_handler;
        
};

#endif /* EXTENDED_REFINE_SCHEDULE */
