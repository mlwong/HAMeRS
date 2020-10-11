/*************************************************************************
 *
 * This file is modified from StandardTagAndInitStrategy.h of the SAMRAI
 * distribution.  For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for error detection.
 *
 ************************************************************************/

#ifndef EXTENDED_TAG_AND_INIT_STRATEGY_HPP
#define EXTENDED_TAG_AND_INIT_STRATEGY_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"

using namespace SAMRAI;

/**
 * Class ExtendedTagAndInitStrategy is an abstract base class
 * that defines a Strategy pattern interface for concrete cell tagging
 * and level initialization operations that are needed by the
 * ExtendedTagAndInitialize class.  This base class insulates
 * that algorithm class from routines for initializing a new level in
 * the hierarchy and for tagging cells to be refined.  Generally, these
 * operations are specific to the problem being solved and the solution
 * methods being employed.  An object of this type is passed into the
 * ExtendedTagAndInitialize constructor.
 *
 * This base class has two pure virtual functions:
 * initializeLevelData(), and resetHierarchyConfiguration() that must
 * be implemented for any regridding method.  The first function sets
 * data on a new level after regridding.  The second function is called
 * at the end of the regridding process and can be used to set communication
 * schedules, for example, which depend on the configuration of the AMR
 * patch hierarchy.  Other routines are virtual here and given default
 * implementations as they are specific to only one type of error estimation
 * method. Value detector functionality requires an implementation of the
 * applyValueDetector() routine. Gradient detector functionality requires an
 * implementation of the applyGradientDetector() routine. Multiresolution
 * detector functionality requires an implementation of the
 * applyMultiresolutionDetector() routine. Integral detector functionality
 * requires an implementation of the applyIntegralDetector() routine. The
 * Richardson extrapolation method requires implementations of the methods:
 * applyRichardsonExtrapolation(), coarsenDataForRichardsonExtrapolation(),
 * getLevelDt(), advanceLevel(), resetTimeDependentData(), and
 * resetDataToPreadvanceState().
 *
 * @see ExtendedTagAndInitialize.
 */

class ExtendedTagAndInitStrategy
{
    public:
        /**
         * Default constructor for
         * ExtendedTagAndInitStrategy.
         */
        ExtendedTagAndInitStrategy();
        
        /**
         * Empty destructor for
         * ExtendedTagAndInitStrategy.
         */
        virtual ~ExtendedTagAndInitStrategy();
        
        /**
         * Determine time increment to advance data on level. The
         * recompute_dt option specifies whether to compute
         * the timestep using the current level data or to return the value
         * stored by the time integrator. The default true setting means
         * the timestep will be computed if no value is supplied.
         *
         * This routine is only when Richardson extrapolation is being used.
         * It is virtual with an empty implementation here (rather than pure
         * virtual) so that users are not required to provide an implementation
         * when the function is not needed.
         */
        virtual double
        getLevelDt(
            const boost::shared_ptr<hier::PatchLevel>& level,
            const double dt_time,
            const bool initial_time);
        
        /**
         * Advance data on all patches on specified patch level from current time
         * (current_time) to new time (new_time).   This routine is called only
         * during time-dependent regridding procedures, such as Richardson
         * extrapolation.  It is virtual with an empty implementation here (rather
         * than pure virtual) so that users are not required to provide an
         * implementation when the function is not needed.  The boolean arguments
         * are used to determine the state of the algorithm and the data when the
         * advance routine is called.  Note that this advance function is also
         * used during normal time integration steps.
         *
         * When this function is called, the level data required to begin the
         * advance must be allocated and be defined appropriately.  Typically,
         * this is equivalent to what is needed to initialize a new level after
         * regridding.  Upon exiting this routine, both current and new data may
         * exist on the level.  This data is needed until level synchronization
         * occurs, in general. Current and new data may be reset by calling
         * the member function resetTimeDependentData().
         *
         * This routine is called from two different points within the Richardson
         * exptrapolation process: to advance a temporary level that is coarser
         * than the hierarchy level on which error estimation is performed, and
         * to advance the hierarchy level itself.  In the first case, the values of
         * the boolean flags are:
         *
         *
         *
         *    - \b  first_step
         *        = true.
         *    - \b  last_step
         *        = true.
         *    - \b  regrid_advance
         *        = true.
         *
         *
         *
         * In the second case, the values of the boolean flags are:
         *
         *
         *
         *    - \b  first_step
         *      (when regridding during time integration sequence)
         *        = true when the level is not coarsest level to synchronize
         *          immediately before the regridding process; else, false.
         *      (when generating initial hierarchy construction)
         *        = true, even though there may be multiple advance steps.
         *    - \b  last_step
         *        = true when the advance is the last in the Richardson
         *          extrapolation step sequence; else false.
         *    - \b  regrid_advance
         *        = true.
         *
         *
         *
         */
        virtual double
        advanceLevel(
            const boost::shared_ptr<hier::PatchLevel>& level,
            const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
            const double current_time,
            const double new_time,
            const bool first_step,
            const bool last_step,
            const bool regrid_advance = false);
        
        /**
         * Reset time-dependent data storage for the specified patch level.
         *
         * This routine only applies when Richardson extrapolation is being used.
         * It is virtual with an empty implementation here (rather than pure
         * virtual) so that users are not required to provide an implementation
         * when the function is not needed.
         */
        virtual void
        resetTimeDependentData(
            const boost::shared_ptr<hier::PatchLevel>& level,
            const double new_time,
            const bool can_be_refined);

        /**
         * Reset data on the patch level by destroying all patch data other
         * than that which is needed to initialize the solution on that level.
         * In other words, this is the data needed to begin a time integration
         * step on the level.
         *
         * This routine is only when Richardson extrapolation is being used.
         * It is virtual with an empty implementation here (rather than pure
         * virtual) so that users are not required to provide an implementation
         * when the function is not needed.
         */
        virtual void
        resetDataToPreadvanceState(
            const boost::shared_ptr<hier::PatchLevel>& level);
        
        /**
         * Initialize data on a new level after it is inserted into an AMR patch
         * hierarchy by the gridding algorithm.  The level number indicates
         * that of the new level.
         *
         * Generally, when data is set, it is interpolated from coarser levels
         * in the hierarchy.  If the old level pointer in the argument list is
         * non-null, then data is copied from the old level to the new level
         * on regions of intersection between those levels before interpolation
         * occurs.   In this case, the level number must match that of the old
         * level.  The specific operations that occur when initializing level
         * data are determined by the particular solution methods in use; i.e.,
         * in the subclass of this abstract base class.
         *
         * The boolean argument initial_time indicates whether the level is
         * being introduced for the first time (i.e., at initialization time),
         * or after some regrid process during the calculation beyond the initial
         * hierarchy construction.  This information is provided since the
         * initialization of the data may be different in each of those
         * circumstances.  The can_be_refined boolean argument indicates whether
         * the level is the finest allowable level in the hierarchy.
         */
        virtual void
        initializeLevelData(
            const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
            const int level_number,
            const double init_data_time,
            const bool can_be_refined,
            const bool initial_time,
            const boost::shared_ptr<hier::PatchLevel>& old_level =
               boost::shared_ptr<hier::PatchLevel>(),
            const bool allocate_data = true) = 0;
        
        /**
         * After hierarchy levels have changed and data has been initialized on
         * the new levels, this routine can be used to reset any information
         * needed by the solution method that is particular to the hierarchy
         * configuration.  For example, the solution procedure may cache
         * communication schedules to amortize the cost of data movement on the
         * AMR patch hierarchy.  This function will be called by the gridding
         * algorithm after the initialization occurs so that the algorithm-specific
         * subclass can reset such things.  Also, if the solution method must
         * make the solution consistent across multiple levels after the hierarchy
         * is changed, this process may be invoked by this routine.  Of course the
         * details of these processes are determined by the particular solution
         * methods in use.
         *
         * The level number arguments indicate the coarsest and finest levels
         * in the current hierarchy configuration that have changed.  It should
         * be assumed that all intermediate levels have changed as well.
         */
        virtual void
        resetHierarchyConfiguration(
            const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
            const int coarsest_level,
            const int finest_level) = 0;
        
        /**
         * Set integer tags to "one" in cells where refinement of the given
         * level should occur according to some user-supplied value criteria.
         * The double time argument is the regrid time.  The integer "tag_index"
         * argument is the patch descriptor index of the cell-centered integer tag
         * array on each patch in the hierarchy.  The boolean argument
         * initial_time indicates whether the level is being subject to refinement
         * at the initial simulation time.  If it is false, then the error
         * estimation process is being invoked at some later time after the AMR
         * hierarchy was initially constructed.  Typically, this information is
         * passed to the user's patch tagging routines since the error
         * estimator or value detector may be different in each case.
         *
         * The boolean uses_gradient_detector_too is true when gradient detector
         * is used in addition to the value detector, and false otherwise.  This
         * argument helps the user to manage multiple regridding criteria.
         *
         * The boolean uses_multiresolution_detector_too is true when
         * multiresolution detector is used in addition to the value detector,
         * and false otherwise.  This argument helps the user to manage multiple
         * regridding criteria.
         *
         * The boolean uses_integral_detector_too is true when integral detector
         * is used in addition to the value detector, and false otherwise.  This
         * argument helps the user to manage multiple regridding criteria.
         * 
         * The boolean uses_richardson_extrapolation_too is true when Richardson
         * extrapolation error estimation is used in addition to the value
         * detector, and false otherwise.  This argument helps the user to
         * manage multiple regridding criteria.
         *
         * This routine is only when value detector is being used.
         * It is virtual with an empty implementation here (rather than pure
         * virtual) so that users are not required to provide an implementation
         * when the function is not needed.
         */
        virtual void
        applyValueDetector(
            const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
            const int level_number,
            const double error_data_time,
            const int tag_index,
            const bool initial_time,
            const bool uses_gradient_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Set integer tags to "one" in cells where refinement of the given
         * level should occur according to some user-supplied gradient criteria.
         * The double time argument is the regrid time.  The integer "tag_index"
         * argument is the patch descriptor index of the cell-centered integer tag
         * array on each patch in the hierarchy.  The boolean argument
         * initial_time indicates whether the level is being subject to refinement
         * at the initial simulation time.  If it is false, then the error
         * estimation process is being invoked at some later time after the AMR
         * hierarchy was initially constructed.  Typically, this information is
         * passed to the user's patch tagging routines since the error
         * estimator or gradient detector may be different in each case.
         *
         * The boolean uses_value_detector_too is true when value detector is
         * used in addition to the gradient detector, and false otherwise.  This
         * argument helps the user to manage multiple regridding criteria.
         *
         * The boolean uses_multiresolution_detector_too is true when
         * multiresolution detector is used in addition to the gradient
         * detector, and false otherwise.  This argument helps the user to
         * manage multiple regridding criteria.
         *
         * The boolean uses_integral_detector_too is true when integral detector
         * is used in addition to the gradient detector, and false otherwise.
         * This argument helps the user to manage multiple regridding criteria.
         * 
         * The boolean uses_richardson_extrapolation_too is true when Richardson
         * extrapolation error estimation is used in addition to the gradient
         * detector, and false otherwise.  This argument helps the user to
         * manage multiple regridding criteria.
         *
         * This routine is only when gradient detector is being used.
         * It is virtual with an empty implementation here (rather than pure
         * virtual) so that users are not required to provide an implementation
         * when the function is not needed.
         */
        virtual void
        applyGradientDetector(
            const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
            const int level_number,
            const double error_data_time,
            const int tag_index,
            const bool initial_time,
            const bool uses_value_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Set integer tags to "one" in cells where refinement of the given
         * level should occur according to some user-supplied multiresolution criteria.
         * The double time argument is the regrid time.  The integer "tag_index"
         * argument is the patch descriptor index of the cell-centered integer tag
         * array on each patch in the hierarchy.  The boolean argument
         * initial_time indicates whether the level is being subject to refinement
         * at the initial simulation time.  If it is false, then the error
         * estimation process is being invoked at some later time after the AMR
         * hierarchy was initially constructed.  Typically, this information is
         * passed to the user's patch tagging routines since the error
         * estimator or multiresolution detector may be different in each case.
         *
         * The boolean uses_value_detector_too is true when value detector is used
         * in addition to the multiresolution detector, and false otherwise.  This
         * argument helps the user to manage multiple regridding criteria.
         *
         * The boolean uses_gradient_detector_too is true when gradient detector is
         * used in addition to the multiresolution detector, and false otherwise.
         * This argument helps the user to manage multiple regridding criteria.
         *
         * The boolean uses_integral_detector_too is true when integral detector
         * is used in addition to the multiresolution detector, and false
         * otherwise.  This argument helps the user to manage multiple regridding
         * criteria.
         * 
         * The boolean uses_richardson_extrapolation_too is true when Richardson
         * extrapolation error estimation is used in addition to the multiresolution
         * detector, and false otherwise.  This argument helps the user to
         * manage multiple regridding criteria.
         * 
         * This routine is only when multiresolution detector is being used.
         * It is virtual with an empty implementation here (rather than pure
         * virtual) so that users are not required to provide an implementation
         * when the function is not needed.
         */
        virtual void
        applyMultiresolutionDetector(
            const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
            const int level_number,
            const double error_data_time,
            const int tag_index,
            const bool initial_time,
            const bool uses_value_detector_too,
            const bool uses_gradient_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Set integer tags to "one" in cells where refinement of the given
         * level should occur according to some user-supplied integral criteria.
         * The double time argument is the regrid time.  The integer "tag_index"
         * argument is the patch descriptor index of the cell-centered integer tag
         * array on each patch in the hierarchy.  The boolean argument
         * initial_time indicates whether the level is being subject to refinement
         * at the initial simulation time.  If it is false, then the error
         * estimation process is being invoked at some later time after the AMR
         * hierarchy was initially constructed.  Typically, this information is
         * passed to the user's patch tagging routines since the error
         * estimator or integral detector may be different in each case.
         *
         * The boolean uses_value_detector_too is true when value detector
         * is used in addition to the integral detector, and false otherwise.
         * This argument helps the user to manage multiple regridding criteria.
         *
         * The boolean uses_multiresolution_detector_too is true when multiresolution
         * detector is used in addition to the integral detector, and false
         * otherwise.  This argument helps the user to manage multiple regridding
         * criteria.
         *
         * The boolean uses_gradient_detector_too is true when gradient detector
         * is used in addition to the integral detector, and false otherwise.
         * This argument helps the user to manage multiple regridding criteria.
         *
         * The boolean uses_multiresolution_detector_too is true when multiresolution
         * detector is used in addition to the integral detector, and false
         * otherwise.  This argument helps the user to manage multiple regridding
         * criteria.
         * 
         * The boolean uses_richardson_extrapolation_too is true when Richardson
         * extrapolation error estimation is used in addition to the integral
         * detector, and false otherwise.  This argument helps the user to
         * manage multiple regridding criteria.
         * 
         * This routine is only when integral detector is being used.
         * It is virtual with an empty implementation here (rather than pure
         * virtual) so that users are not required to provide an implementation
         * when the function is not needed.
         */
        virtual void
        applyIntegralDetector(
            const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
            const int level_number,
            const double error_data_time,
            const int tag_index,
            const bool initial_time,
            const bool uses_value_detector_too,
            const bool uses_gradient_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Set integer tags to "one" in cells where refinement of the given
         * level should occur according to some user-supplied Richardson
         * extrapolation criteria.  The "error_data_time" argument is the
         * regrid time.  The "deltat" argument is the time increment to advance
         * the solution on the level to be refined.  Note that that level is
         * finer than the level in the argument list, in general.  The
         * ratio between the argument level and the actual hierarchy level
         * is given by the integer "coarsen ratio".
         *
         * The integer "tag_index" argument is the patch descriptor index of
         * the cell-centered integer tag array on each patch in the hierarchy.
         *
         * The boolean argument initial_time indicates whether the level is being
         * subject to refinement at the initial simulation time.  If it is false,
         * then the error estimation process is being invoked at some later time
         * after the AMR hierarchy was initially constructed.  Typically, this
         * information is passed to the user's patch tagging routines since the
         * application of the Richardson extrapolation process may be different
         * in each case.
         *
         * The boolean uses_value_detector_too is true when a value detector
         * is used in addition to Richardson extrapolation, and false otherwise.
         * This argument helps the user to manage multiple regridding criteria.
         *
         * The boolean uses_gradient_detector_too is true when a gradient detector
         * is used in addition to Richardson extrapolation, and false otherwise.
         * This argument helps the user to manage multiple regridding criteria.
         *
         * The boolean uses_multiresolution_detector_too is true when a
         * multiresolution detector is used in addition to the Richardson
         * extrapolation, and false otherwise.  This argument helps the user to
         * manage multiple regridding criteria.
         *
         * The boolean uses_integral_detector_too is true when integral detector
         * is used in addition to the Richardson extrapolation, and false
         * otherwise.  This argument helps the user to manage multiple regridding
         * criteria.
         * 
         * This routine is only when Richardson extrapolation is being used.
         * It is virtual with an empty implementation here (rather than pure
         * virtual) so that users are not required to provide an implementation
         * when the function is not needed.
         */
        virtual void
        applyRichardsonExtrapolation(
            const boost::shared_ptr<hier::PatchLevel>& level,
            const double error_data_time,
            const int tag_index,
            const double deltat,
            const int error_coarsen_ratio,
            const bool initial_time,
            const bool uses_value_detector_too,
            const bool uses_gradient_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_integral_detector_too);
        
        /**
         * Coarsen solution data from level to coarse_level for Richardson
         * extrapolation.  Note that this routine will be called twice during
         * the Richardson extrapolation error estimation process, once to set
         * data on the coarser level and once to coarsen data from after
         * advancing the fine level.  The init_coarse_level boolean argument
         * indicates whether data is set on the coarse level by coarsening the
         * "old" time level solution or by coarsening the "new" solution on the
         * fine level (i.e., after it has been advanced).
         *
         * This routine is only when Richardson extrapolation is being used.
         * It is virtual with an empty implementation here (rather than pure
         * virtual) so that users are not required to provide an implementation
         * when the function is not needed.
         */
        virtual void
        coarsenDataForRichardsonExtrapolation(
           const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
           const int level_number,
           const boost::shared_ptr<hier::PatchLevel>& coarser_level,
           const double coarsen_data_time,
           const bool before_advance);
        
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
         * This is a virtual method with a no-op implementation provided, so that
         * users who do not need this processing step need not implement anything.
         *
         * @param hierarchy The PatchHierarchy being modified.
         * @param level_number The number of the PatchLevel in hierarchy being
         *                     added or regridded.
         * @param new_box_level BoxLevel containing the boxes for the new level
         *
         */
        virtual void
        processHierarchyBeforeAddingNewLevel(
           const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
           const int level_number,
           const boost::shared_ptr<hier::BoxLevel>& new_box_level);
        
        /**
         * In some cases user code may wish to process a PatchLevel before it is
         * removed from the hierarchy.  For example, data may exist only on a given
         * PatchLevel such as the finest level.  If that level were to be removed
         * before this data is moved off of it then the data will be lost.  This
         * method is a user defined callback used by GriddingAlgorithm when a
         * PatchLevel is to be removed.  The callback performs any user actions on
         * the level about to be removed.  It is implemented by classes derived from
         * ExtendedTagAndInitStrategy.
         *
         * This is a virtual method with a no-op implementation provided, so that
         * users who do not need this processing step need not implement anything.
         *
         * @param hierarchy The PatchHierarchy being modified.
         * @param level_number The number of the PatchLevel in hierarchy about to be
         *                     removed.
         * @param old_level The level in hierarchy about to be removed.
         *
         * @see mesh::GriddingAlgorithm
         * @see ExtendedTagAndInitStrategy
         */
        virtual void
        processLevelBeforeRemoval(
            const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
            const int level_number,
            const boost::shared_ptr<hier::PatchLevel>& old_level =
                boost::shared_ptr<hier::PatchLevel>());
        
    private:
        
};

#endif /* EXTENDED_TAG_AND_INIT_STRATEGY_HPP */

