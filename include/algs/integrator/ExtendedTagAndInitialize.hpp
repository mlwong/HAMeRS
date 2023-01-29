/*************************************************************************
 *
 * This file is modified from StandardTagAndInitialize.h of the SAMRAI
 * distribution.  For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Routines for performing cell-tagging and initializing
 *                a new level.
 *
 ************************************************************************/

#ifndef EXTENDED_TAG_AND_INITIALIZE_HPP
#define EXTENDED_TAG_AND_INITIALIZE_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "algs/integrator/ExtendedTagAndInitializeConnectorWidthRequestor.hpp"
#include "algs/integrator/ExtendedTagAndInitStrategy.hpp"

#include "SAMRAI/mesh/TagAndInitializeStrategy.h"

#include <vector>

using namespace SAMRAI;

/*!
 * Class ExtendedTagAndInitialize defines an implementation
 * for level initialization and cell tagging routines needed by
 * the GriddingAlgorithm class.  This class is derived from
 * the abstract base class TagAndInitializeStrategy.  It invokes
 * problem-specific level initialization routines after AMR patch
 * hierarchy levels change and routines for tagging cells for refinement
 * using one (or more) of the following methods:
 *
 *   - Value Detection
 *   - Gradient Detection
 *   - Multiresolution Detection
 *   - Integral Detection
 *   - Richardson Extrapolation
 *   - Explicitly defined refine boxes
 *
 * Tagging methods may be activated at specific cycles and/or times.
 * It is possible to use combinations of the tagging methods (e.g., use
 * value detection, gradient detection, multiresolution detection, Richardson
 * extrapolation, and static refine boxes at the same cycle/time).  The order
 * in which they are executed is fixed (Richardson extrapolation first,
 * multiresolution detection second, gradient detection third, value detection
 * fourth and refine boxes last). An input entry for this class is optional.
 * If none is provided, the class will, by default, not use any criteria
 * to tag cells for refinement and issue a warning.
 *
 * <b> Input Parameters </b>
 *
 * <b> Definitions: </b> <br>
 * The general input syntax is as follows:
 *
 *    - \b    at_0
 *      input section describing a set of tagging methods to apply to a
 *      cycle or time
 *         - \b cycle = integer cycle at which the set of tagging methods is
 *                      activated
 *                      one of cycle or time must be supplied
 *         - \b time = double time at which the set of tagging methods is
 *                     activated
 *                     one of cycle or time must be supplied
 *         - \b tag_0
 *           first tagging method in this set of tagging methods
 *              - \b tagging_method = one of RICHARDSON_EXTRAPOLATION,
 *                                    MULTIRESOLUTION_DETECTOR, GRADIENT_DETECTOR
 *                                    VALUE_DETECTOR, REFINE_REGIONS, REFINE_BOXES, NONE
 *              - \b level_m
 *                required if tagging_method is REFINE_BOXES, the static boxes
 *                for the mth level
 *                   - \b block_m
 *                     required description of refine boxes in the mth block
 *                     of the level just specified
 *                        - \b boxes = required box array specifying refine
 *                                     boxes
 *                   - \b . . .
 *                   - \b block_n
 *              - \b . . .
 *              - \b level_n
 *         - \b . . .
 *         - \b tag_n
 *    - \b . . .
 *    - \b    at_n
 *
 *       If both @b time or @b cycle entries are supplied, at any point
 *       where both a time and a cycle entry are active, the time entry
 *       takes precedence.
 *
 *       If the problem is single block one may omit the "block_0" database
 *       enclosing the boxes and simply include the boxes directly inside the
 *       "level_m" database.
 *
 *       It is possible to use a "shortcut" input syntax for extremely
 *       simple tagging criteria.  If you only want RICHARDSON_EXTRAPOLATION,
 *       MULTIRESOLUTION_DETECTOR, GRADIENT_DETECTOR, VALUE_DETECTOR, REFINE_REGIONS on for
 *       the entire simulation then an input of the following form may be used:
 *
 * @code
 *    tagging_method = RICHARDSON_EXTRAPOLATION
 * @endcode
 *
 * A sample input file entry for a multi-block problem might look like:
 *
 * @code
 *    at_0 {
 *       cycle = 0
 *       tag_0 {
 *          tagging_method = REFINE_BOXES
 *          level_0 {
 *             block_1 {
 *                boxes = [(5,5),(9,9)],[(12,15),(18,19)]
 *             }
 *          }
 *          level_1 {
 *             block_1 {
 *                boxes = [(25,30),(29,35)]
 *             }
 *          }
 *          level_2 {
 *             block_1 {
 *                boxes = [(60,70),(70,80)]
 *             }
 *          }
 *       }
 *    }
 *
 *    at_1 {
 *       cycle = 10
 *       tag_0 {
 *          tagging_method = REFINE_BOXES
 *          level_0 {
 *             block_2 {
 *                boxes = [(7,7),(11,11)],[(14,17),(20,21)]
 *             }
 *          }
 *       }
 *    }
 *
 *    at_2 {
 *       time = 0.05
 *       tag_0 {
 *          tagging_method = REFINE_BOXES
 *          level_1 {
 *             block_0 {
 *                boxes = [(30,35),(34,40)]
 *             }
 *          }
 *       }
 *    }
 *
 *    at_3 {
 *       time = 0.10
 *       tag_0 {
 *          tagging_method = REFINE_BOXES
 *          level_1 {
 *             block_1 {
 *                boxes = [(35,40),(39,45)]
 *             }
 *          }
 *       }
 *    }
 * @endcode
 *
 * This class supplies the routines for tagging cells
 * and invoking problem-specific level initialization routines after AMR
 * patch hierarchy levels change.  A number of problem-specific operations
 * are performed in the ExtendedTagAndInitStrategy
 * data member, for which methods are specified in a derived subclass.
 *
 * @see mesh::TagAndInitializeStrategy
 * @see mesh::GriddingAlgorithm
 * @see ExtendedTagAndInitStrategy
 */

class ExtendedTagAndInitialize:
    public mesh::TagAndInitializeStrategy
{
    public:
        /*!
         * Constructor for ExtendedTagAndInitialize which
         * may read inputs from the provided input_db.  If no input
         * database is provided, the class interprets that no tagging
         * is desired so no cell-tagging will be performed.
         *
         * @pre !object_name.empty()
         */
        ExtendedTagAndInitialize(
            const std::string& object_name,
            ExtendedTagAndInitStrategy* tag_strategy,
            const HAMERS_SHARED_PTR<tbox::Database>& input_db =
                HAMERS_SHARED_PTR<tbox::Database>());
        
        /*!
         * Virtual destructor for ExtendedTagAndInitialize.
         */
        virtual ~ExtendedTagAndInitialize();
        
        /*!
         * Specifies whether the chosen method advances the solution data
         * in the regridding process at any cycle or time (Richardson
         * extrapolation does, the others will not).
         */
        bool
        everUsesTimeIntegration() const;
        
        /*!
         * Specifies whether the chosen method advances the solution data
         * in the regridding process at the supplied cycle or time (Richardson
         * extrapolation does, the others will not).
         */
        bool
        usesTimeIntegration(
            int cycle,
            double time);
        
        /*!
         * Returns true if Richardson extrapolation is used at any cycle or
         * time.
         */
        bool
        everUsesRichardsonExtrapolation() const;
        
        /*!
         * Returns true if Richardson extrapolation is used at the supplied cycle
         * and time.
         *
         * @pre !d_use_cycle_criteria || !d_use_time_criteria
         */
        bool
        usesRichardsonExtrapolation(
            int cycle,
            double time);
        
        /*!
         * Returns true if integral detector is used at any cycle or time.
         */
        bool
        everUsesIntegralDetector() const;
        
        /*!
         * Returns true if integral detector is used at the supplied cycle and time.
         *
         * @pre !d_use_cycle_criteria || !d_use_time_criteria
         */
        bool
        usesIntegralDetector(
            int cycle,
            double time);
        
        /*!
         * Returns true if multiresolution detector is used at any cycle or time.
         */
        bool
        everUsesMultiresolutionDetector() const;
        
        /*!
         * Returns true if multiresolution detector is used at the supplied cycle and time.
         *
         * @pre !d_use_cycle_criteria || !d_use_time_criteria
         */
        bool
        usesMultiresolutionDetector(
            int cycle,
            double time);
        
        /*!
         * Returns true if gradient detector is used at any cycle or time.
         */
        bool
        everUsesGradientDetector() const;
        
        /*!
         * Returns true if gradient detector is used at the supplied cycle and time.
         *
         * @pre !d_use_cycle_criteria || !d_use_time_criteria
         */
        bool
        usesGradientDetector(
            int cycle,
            double time);
        
        /*!
         * Returns true if value detector is used at any cycle or time.
         */
        bool
        everUsesValueDetector() const;
        
        /*!
         * Returns true if value detector is used at the supplied cycle and time.
         *
         * @pre !d_use_cycle_criteria || !d_use_time_criteria
         */
        bool
        usesValueDetector(
            int cycle,
            double time);
        
        /*!
         * Returns true if user supplied refine regions are used at any cycle or time.
         */
        bool
        everUsesRefineRegions() const;
        
        /*!
         * Returns true if user supplied refine regions are used at the supplied cycle and time.
         *
         * @pre !d_use_cycle_criteria || !d_use_time_criteria
         */
        bool
        usesRefineRegions(
            int cycle,
            double time);
        
        /*!
         * Returns true if user supplied refine boxes are used at any cycle or time.
         */
        bool
        everUsesRefineBoxes() const;
        
        /*!
         * Returns true if user supplied refine boxes are used at the supplied cycle
         * and time.
         *
         * @pre !d_use_cycle_criteria || !d_use_time_criteria
         */
        bool
        usesRefineBoxes(
           int cycle,
           double time);
        
        /*!
         * Return coarsen ratio used for applying cell tagging. An error
         * coarsen ratio other than 2 or 3 will throw an error.
         */
        int
        getErrorCoarsenRatio() const;
        
        /*!
         * Some restrictions may be placed on the coarsen ratio used for
         * cell tagging.  Check these here.
         */
        void
        checkCoarsenRatios(
            const std::vector<hier::IntVector>& ratio_to_coarser);
        
        /*!
         * Pass the request to initialize the data on a new level in the
         * hierarchy to the ExtendedTagAndInitStrategy data member. Required
         * arguments specify the grid hierarchy, level number being initialized,
         * simulation time at which the data is initialized, whether the level
         * can be refined, and whether it is the initial time.  Optional arguments
         * include an old level, from which data may be used to initialize this
         * level, and a flag that indicates whether data on the initialized level
         * must first be allocated.  For more information on the operations that
         * must be performed, see the
         * TagAndInitializeStrategy::initializeLevelData() method.
         *
         * @pre hierarchy
         * @pre (level_number >= 0) &&
         *      (level_number <= hierarchy->getFinestLevelNumber())
         * @pre !old_level || (level_number == old_level->getLevelNumber())
         * @pre hierarchy->getPatchLevel(level_number)
         */
        void
        initializeLevelData(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
            const int level_number,
            const double init_data_time,
            const bool can_be_refined,
            const bool initial_time,
            const HAMERS_SHARED_PTR<hier::PatchLevel>& old_level =
               HAMERS_SHARED_PTR<hier::PatchLevel>(),
            const bool allocate_data = true);
        
        /*!
         * Pass the request to reset information that depends on the hierarchy
         * configuration to the ExtendedTagAndInitStrategy data member.
         * For more information on the operations that must be performed, see
         * the TagAndInitializeStrategy::resetHierarchyConfiguration()
         * method.
         *
         * @pre hierarchy
         * @pre (coarsest_level >= 0) && (coarsest_level <= finest_level) &&
         *      (finest_level <= hierarchy->getFinestLevelNumber())
         */
        void
        resetHierarchyConfiguration(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
            const int coarsest_level,
            const int finest_level);
        
        /*!
         * Certain cases may require pre-processing of error estimation data
         * before tagging cells, which is handled by this method.  For more
         * information on the operations that must be performed, see the
         * TagAndInitializeStrategy::preprocessErrorEstimation()
         * method.
         *
         * @pre hierarchy
         * @pre (level_number >= 0) &&
         *      (level_number <= hierarchy->getFinestLevelNumber())
         * @pre hierarchy->getPatchLevel(level_number)
         */
        void
        preprocessErrorEstimation(
           const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
           const int level_number,
           const int cycle,
           const double regrid_time,
           const double regrid_start_time,
           const bool initial_time);
        
        /*!
         * Pass the request to set tags on the given level where refinement of
         * that level should occur.  Gradient detection, Richardson extrapolation,
         * and tagging on static refine boxes is performed here.
         *
         * For more information on the operations that must be performed, see the
         * TagAndInitializeStrategy::tagCellsForRefinement() routine.
         *
         * @pre level
         * @pre (level_number >= 0) &&
         *      (level_number <= level->getFinestLevelNumber())
         * @pre level->getPatchLevel(level_number)
         * @pre tag_index >= 0
         */
        void
        tagCellsForRefinement(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& level,
            const int level_number,
            const int regrid_cycle,
            const double regrid_time,
            const int tag_index,
            const bool initial_time,
            const bool coarsest_sync_level,
            const bool can_be_refined,
            const double regrid_start_time = 0);
        
        /*!
         * Return true if boxes for coarsest hierarchy level are not appropriate
         * for gridding strategy.  Otherwise, return false.  If false is returned,
         * it is useful to provide a detailed explanatory message describing the
         * problems with the boxes.
         *
         * @pre !boxes.empty()
         */
        bool
        coarsestLevelBoxesOK(
            const hier::BoxContainer& boxes) const;
        
        /*!
         * Return whether refinement is being performed using ONLY
         * user-supplied refine boxes.  If any method is used that invokes
         * tagging, this will return false.
         */
        bool
        refineUserBoxInputOnly(
            int cycle,
            double time);
        
        const ExtendedTagAndInitializeConnectorWidthRequestor&
        getConnectorWidthRequestor() const
        {
            return d_staicwri;
        }
        
        /*!
         * Return user supplied set of refine boxes for specified level number
         * and time.  The boolean return value specifies whether the boxes
         * have been reset from the last time this method was called.  If they
         * have been reset, it returns true.  If they are unchanged, it returns
         * false.
         *
         * @pre refine_boxes.empty()
         * @pre level_num >= 0
         * @pre time >= 0.0
         * @pre !d_use_cycle_criteria || !d_use_time_criteria
         */
        bool
        getUserSuppliedRefineBoxes(
            hier::BoxContainer& refine_boxes,
            const int level_number,
            const int cycle,
            const double time);
        
        /*!
         * Reset the static refine boxes for the specified level number in the
         * hierarchy.  The level number must be greater than or equal to zero.
         *
         * @pre level_num >= 0
         */
        void
        resetRefineBoxes(
            const hier::BoxContainer& refine_boxes,
            const int level_number);
        
        /*!
         * Turn on refine boxes criteria at the specified time programmatically.
         *
         * @param time Time to turn refine boxes criteria on.
         */
        void
        turnOnRefineBoxes(
            double time);
        
        /*!
         * Turn off refine boxes criteria at the specified time programmatically.
         *
         * @param time Time to turn refine boxes criteria off.
         */
        void
        turnOffRefineBoxes(
            double time);
        
        /*!
         * Turn on refine regions criteria at the specified time
         * programmatically.
         *
         * @param time Time to turn refine regions criteria on.
         *
         * @pre d_tag_strategy
         */
        void
        turnOnRefineRegions(
            double time);
        
        /*!
         * Turn off refine regions criteria at the specified time
         * programmatically.
         *
         * @param time Time to turn refine regions criteria off.
         */
        void
        turnOffRefineRegions(
            double time);
        
        /*!
         * Turn on value detector criteria at the specified time
         * programmatically.
         *
         * @param time Time to turn value detector criteria on.
         *
         * @pre d_tag_strategy
         */
        void
        turnOnValueDetector(
            double time);
        
        /*!
         * Turn off value detector criteria at the specified time
         * programmatically.
         *
         * @param time Time to turn value detector criteria off.
         */
        void
        turnOffValueDetector(
            double time);
        
        /*!
         * Turn on gradient detector criteria at the specified time
         * programmatically.
         *
         * @param time Time to turn gradient detector criteria on.
         *
         * @pre d_tag_strategy
         */
        void
        turnOnGradientDetector(
            double time);
        
        /*!
         * Turn off gradient detector criteria at the specified time
         * programmatically.
         *
         * @param time Time to turn gradient detector criteria off.
         */
        void
        turnOffGradientDetector(
            double time);
        
        /*!
         * Turn on multiresolution detector criteria at the specified time
         * programmatically.
         *
         * @param time Time to turn multiresolution detector criteria on.
         *
         * @pre d_tag_strategy
         */
        void
        turnOnMultiresolutionDetector(
            double time);
        
        /*!
         * Turn off multiresolution detector criteria at the specified time
         * programmatically.
         *
         * @param time Time to turn multiresolution detector criteria off.
         */
        void
        turnOffMultiresolutionDetector(
            double time);
        
        /*!
         * Turn on integral detector criteria at the specified time
         * programmatically.
         *
         * @param time Time to turn integral detector criteria on.
         *
         * @pre d_tag_strategy
         */
        void
        turnOnIntegralDetector(
            double time);
        
        /*!
         * Turn off integral detector criteria at the specified time
         * programmatically.
         *
         * @param time Time to turn integral detector criteria off.
         */
        void
        turnOffIntegralDetector(
            double time);
        
        /*!
         * Turn on Richardson extrapolation criteria at the specified time
         * programmatically.
         *
         * @param time Time to turn Richardson extrapolation on.
         *
         * @pre d_tag_strategy
         */
        void
        turnOnRichardsonExtrapolation(
            double time);
        
        /*!
         * Turn off Richardson extrapolation criteria at the specified time
         * programmatically.
         *
         * @param time Time to turn Richardson extrapolation off.
         */
        void
        turnOffRichardsonExtrapolation(
            double time);
        
        /*!
         * @brief Process a hierarchy before swapping old and new levels during
         * regrid.
         *
         * During regrid, if user code needs to do any application-specific
         * operations on the PatchHierarchy before a new level is added or
         * an old level is swapped for a new level, this method provides a callback
         * for the user to define such operations.  The PatchHierarchy is provided
         * in its state with the old level, if it exists, still in place, while the
         * new BoxLevel is also provided so that the user code can know the boxes
         * that will make up the new level.
         *
         * @param hierarchy The PatchHierarchy being modified.
         * @param level_number The number of the PatchLevel in hierarchy being
         *                     added or regridded.
         * @param new_box_level BoxLevel containing the boxes for the new level
         *
         */
        virtual void
        processHierarchyBeforeAddingNewLevel(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
            const int level_number,
            const HAMERS_SHARED_PTR<hier::BoxLevel>& new_box_level);
        
        /*!
         * In some cases user code may wish to process a PatchLevel before it is
         * removed from the hierarchy.  For example, data may exist only on a given
         * PatchLevel such as the finest level.  If that level were to be removed
         * before this data is moved off of it then the data will be lost.  This
         * method is a user defined callback used by GriddingAlgorithm when a
         * PatchLevel is to be removed.  The callback performs any user actions on
         * the level about to be removed.  It is implemented by classes derived from
         * ExtendedTagAndInitStrategy.
         *
         * @param hierarchy The PatchHierarchy being modified.
         * @param level_number The number of the PatchLevel in hierarchy about to be
         *                     removed.
         * @param old_level The level in hierarchy about to be removed.
         *
         * @see mesh::GriddingAlgorithm
         * @see ExtendedTagAndInitStrategy
         */
        void
        processLevelBeforeRemoval(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
            const int level_number,
            const HAMERS_SHARED_PTR<hier::PatchLevel>& old_level =
                HAMERS_SHARED_PTR<hier::PatchLevel>());
        
    private:
        /*
         * Apply preprocessing for Richardson extrapolation.
         */
        void
        preprocessRichardsonExtrapolation(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
            const int level_number,
            const double regrid_time,
            const double regrid_start_time,
            const bool initial_time);
        
        /*
         * Apply Richardson extrapolation algorithm.
         */
        void
        tagCellsUsingRichardsonExtrapolation(
           const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
           const int level_number,
           const int regrid_cycle,
           const double regrid_time,
           const double regrid_start_time,
           const int tag_index,
           const bool initial_time,
           const bool coarsest_sync_level,
           const bool can_be_refined);
        
        /*
         * Set the tagging criteria for the supplied time and cycle.  If there
         * are both cycle and time criteria for this cycle/time combination the
         * time criteria has precedence.
         */
        void
        setCurrentTaggingCriteria(
            int cycle,
            double time);
        
        /*
         * Read input values, indicated above, from given database.
         */
        void
        getFromInput(
           const HAMERS_SHARED_PTR<tbox::Database>& input_db);
        
        /*
         * Concrete object that supplies problem-specific initialization
         * and regridding operations.
         */
        ExtendedTagAndInitStrategy* d_tag_strategy;
        
        /*
         * The error_coarsen_ratio used for all levels in the hierarchy.
         * If Richardson extrapolation is not used, the error coarsen ratio
         * is 1.  If Richardson extrapolation is used, the error coarsen ratio
         * is set in the method coarsestLevelBoxesOK().
         */
        int d_error_coarsen_ratio;
        
        /*
         * tbox::Array of patch levels containing coarsened versions of the patch
         * levels, for use with Richardson extrapolation.
         */
        std::vector<HAMERS_SHARED_PTR<hier::PatchLevel> > d_rich_extrap_coarsened_levels;
        
        ExtendedTagAndInitializeConnectorWidthRequestor d_staicwri;
        
        /*
         * Arrays to hold boxes that are specifically reset by the user (via the
         * resetRefineBoxes() method).  The boolean array specifies which levels
         * have been reset while the box array specifies the new set of refine
         * boxes for the level.
         */
        std::vector<bool> d_refine_boxes_reset;
        std::vector<hier::BoxContainer> d_reset_refine_boxes;
        
        struct TagCriteria
        {
            std::string d_tagging_method;
            std::map<int, hier::BoxContainer> d_level_refine_boxes;
        };
        
        struct CycleTagCriteria
        {
            int d_cycle;
            std::vector<TagCriteria> d_tag_criteria;
        };
        
        struct TimeTagCriteria
        {
            double d_time;
            std::vector<TagCriteria> d_tag_criteria;
        };
        
        struct cycle_tag_criteria_less
        {
            bool
            operator () (const CycleTagCriteria& c1, const CycleTagCriteria& c2) const
            {
                return c1.d_cycle < c2.d_cycle;
            }
        };
        
        struct time_tag_criteria_less
        {
            bool
            operator () (const TimeTagCriteria& c1, const TimeTagCriteria& c2) const
            {
                return c1.d_time < c2.d_time;
            }
        };
        
        /*
         * User defined tagging criteria to be applied at specific cycles ordered
         * by the application cycle.
         */
        std::set<CycleTagCriteria, cycle_tag_criteria_less> d_cycle_criteria;
        
        /*
         * User defined tagging criteria to be applied at specific times ordered
         * by the application time.
         */
        std::set<TimeTagCriteria, time_tag_criteria_less> d_time_criteria;
        
        /*
         * Iterator pointing to the cycle tagging criteria in use if
         * d_use_cycle_criteria is true.
         */
        std::set<CycleTagCriteria, cycle_tag_criteria_less>::iterator
            d_cur_cycle_criteria;
        
        /*
         * Iterator pointing to the time tagging criteria in use if
         * d_use_time_criteria is true.
         */
        std::set<TimeTagCriteria, time_tag_criteria_less>::iterator
            d_cur_time_criteria;
        
        /*
         * Flag indicating if a cycle tagging criteria is now being used.  At most
         * one of d_use_cycle_criteria and d_use_time_criteria can be true.
         */
        bool d_use_cycle_criteria;
        
        /*
         * Flag indicating if a time tagging criteria is now being used.  At most
         * one of d_use_cycle_criteria and d_use_time_criteria can be true.
         */
        bool d_use_time_criteria;
        
        /*
         * Flag indicating if any tagging criteria is RICHARDSON_EXTRAPOLATION.
         */
        bool d_ever_uses_richardson_extrapolation;
        
        /*
         * Flag indicating if any tagging criteria is INTEGRAL_DETECTOR.
         */
        bool d_ever_uses_integral_detector;
        
        /*
         * Flag indicating if any tagging criteria is MULTIRESOLUTION_DETECTOR.
         */
        bool d_ever_uses_multiresolution_detector;
        
        /*
         * Flag indicating if any tagging criteria is GRADIENT_DETECTOR.
         */
        bool d_ever_uses_gradient_detector;
        
        /*
         * Flag indicating if any tagging criteria is VALUE_DETECTOR.
         */
        bool d_ever_uses_value_detector;
        
        /*
         * Flag indicating if any tagging criteria is REFINE_REGIONS.
         */
        bool d_ever_uses_refine_regions;
        
        /*
         * Flag indicating if any tagging criteria is REFINE_BOXES.
         */
        bool d_ever_uses_refine_boxes;
        
        /*
         * Flag indicating if the tagged boxes have changed due to a change in the
         * tagging criteria.
         */
        bool d_boxes_changed;
        
        /*
         * New tagging criteria are set only when the cycle advances.  This tracks
         * the last cycle so the code can tell when the cycle has advanced.
         */
        int d_old_cycle;
        
};

#endif /* EXTENDED_TAG_AND_INITIALIZE_HPP */
