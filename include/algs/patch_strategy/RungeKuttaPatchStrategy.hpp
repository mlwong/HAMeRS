/*************************************************************************
 *
 * This file is modified from HyperbolicPatchStrategy.h of the SAMRAI
 * distribution. For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Interface to patch routines for hyperbolic Runge-Kutta
 *                integration scheme.
 *
 ************************************************************************/


#ifndef RUNGE_KUTTA_PATCH_STRATEGY_HPP
#define RUNGE_KUTTA_PATCH_STRATEGY_HPP

#include "HAMeRS_config.hpp"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include "boost/shared_ptr.hpp"

class RungeKuttaLevelIntegrator;

/**
 * Class RungeKuttaPatchStrategy is an abstract base class defining the interface between an
 * RungeKuttaLevelIntegrator object and operations applied to a single patch in a structured AMR
 * hierarchy. The operations include patch initialization, dt calculation, hyperbolic flux computation,
 * conservative differencing, and error estimation. This class is derived from the
 * xfer::RefinePatchStrategy and xfer::CoarsenPatchStrategy abstract base classes. These base classes
 * provide the interface for user-defined interlevel data refining and coarsening operations and the
 * specification of physical boundary conditions. The functions setPhysicalBoundaryConditions(), and
 * pre/postprocessRefine() are overloaded from the class xfer::RefinePatchStrategy. The operations
 * pre/postprocessCoarsen() are overloaded from xfer::CoarsenPatchStrategy. The pre/postprocessCoarsen/
 * Refine() operations are given empty implementations here so that the user does not need to proovide
 * them if the operations are not needed.
 *
 * It is important to recognize that for the concrete patch strategy subclass and the
 * RungeKuttaLevelIntegrator to work together, the concrete strategy must know which patch data to
 * operate on. The patch data storage is manipulated by the level integrator. The set/clearDataContext()
 * methods allow the integrator to inform the patch strategy of the correct data context. The concrete
 * patch strategy subclass can access the appropriate context via the getDataContext() method.
 *
 * @see RungeKuttaLevelIntegrator
 * @see xfer::RefinePatchStrategy
 * @see xfer::CoarsenPatchStrategy
 */

using namespace SAMRAI;

class RungeKuttaPatchStrategy:
   public xfer::RefinePatchStrategy,
   public xfer::CoarsenPatchStrategy
{
    public:
        /**
         * Blank constructor for RungeKuttaPatchStrategy.
         */
        RungeKuttaPatchStrategy();
        
        /**
         * Virtual destructor for RungeKuttaPatchStrategy.
         */
        virtual ~RungeKuttaPatchStrategy();
        
        /**
         * Register specific variables needed in the numerical routines with the Runge-Kutta level
         * integrator using the registerVariable() function in that class. The integrator manipulates
         * storage for the data and this registration defines the way in which data for each quantity
         * will be manipulated on the patches. Typically, the derived data quantities for plotting
         * are registered with a visualization data writer in this routine as well, since the Runge-
         * Kutta level integrator provides the variable context for plotting (i.e., which data is
         * available when a plot file is generated). The integrator pointer cannot be null in most
         * cases.
         *
         * The gridding algorithm pointer is provided so that patch data objects may be registered
         * with the load balancer object (owned by the gridding algorithm) for non-uniform load
         * balancing, if needed.
         */
        virtual void
        registerModelVariables(
            RungeKuttaLevelIntegrator* integrator) = 0;
        
        /**
         * Set up parameters in the load balancer object (owned by the gridding algorithm) if needed.
         * This function is called immediately after the registerModelVariables() function is called.
         * The Runge-Kutta level integrator pointer is provided so that the integrator can be used
         * to manage data for the load balancer if needed (e.g., when using non-uniform load balancing).
         *
         * Note that this function is not pure virtual. It is given a dummy implementation here so
         * that users may ignore it when inheriting from this class.
         */
        virtual void
        setupLoadBalancer(
            RungeKuttaLevelIntegrator* integrator,
            mesh::GriddingAlgorithm* gridding_algorithm);
        
        /**
         * Set the initial data on a patch interior only. Note that no ghost cells need to be set in
         * this routine regardless of whether the patch data corresponding to the context requires
         * ghost cells. The data_time is the simulation time when the routine is called. The boolean
         * initial_time is true if the routine is called at the initial time when the hierarchy is
         * initially constructed, otherwise it is false.
         */
        virtual void
        initializeDataOnPatch(
            hier::Patch& patch,
            const double data_time,
            const bool initial_time) = 0;
        
        /**
         * Compute the stable time increment for a patch on the level with the given number. The
         * boolean flag initial_time is true if the routine is called at the initial simulation
         * time; otherwise it is false. The double argument dt_time is the simulation time.
         */
        virtual double
        computeStableDtOnPatch(
            hier::Patch& patch,
            const bool initial_time,
            const double dt_time) = 0;
        
        /**
         * Compute TIME INTEGRALS of fluxes to be used in finite difference for patch integration.
         * That is, it is assumed that this numerical routine will compute the fluxes corresponding
         * to the cell faces multiplied by the time increment. Typically, the numerical flux is the
         * normal flux at the cell face. The flux integrals will be used in the computing intermediate
         * fluxes that updates the conserved quantities in advanceSingleStepOnPatch().
         *
         * If the equations have source terms, the TIME INTEGRALS of the sources are also computed.
         *
         * Note that the numerical routines in this method generally require ghost cells. Ghost cells
         * data is filled before this routine is called.
         */
        virtual void
        computeFluxesAndSourcesOnPatch(
            hier::Patch& patch,
            const double time,
            const double dt,
            const int RK_step_number,
            const boost::shared_ptr<hier::VariableContext>& data_context = boost::shared_ptr<hier::VariableContext>()) = 0;      
        
        /**
         * Advance a single Runge-Kutta step using previous intermediate solutions and fluxes computed
         * in computeFluxesAndSourcesOnPatch() routine. Note that the
         * computeFluxesAndSourcesOnPatch() routine computes TIME INTEGRALs of the numerical
         * fluxes and sources (e.g., they have been multiplied by dt) so the time-stepping routine
         * should be consistent with that.
         *
         * @param patch                 patch that RK step is being applied
         * @param time                  simulation time when the routine is called
         * @param dt                    timestep
         * @param alpha                 vector of coefficients applied to the intermediate solution
         *                              in the RK step
         * @param beta                  vector of coefficients applied to the flux
         * @param gamma                 vector of coefficients for accumulation of flux
         * @param intermediate_context  vector of data context for intermediate solution/ flux at
         *                              different stages
         */
        virtual void
        advanceSingleStepOnPatch(
            hier::Patch& patch,
            const double time,
            const double dt,
            const std::vector<double>& alpha,
            const std::vector<double>& beta,
            const std::vector<double>& gamma,
            const std::vector<boost::shared_ptr<hier::VariableContext> >& intermediate_context)
                = 0;
        
        /**
         * Correct the fluxes at the coarse-fine boundaries during a flux synchronization step. Note
         * that the computeFluxesAndSourcesOnPatch() routine computes TIME INTEGRALs of the numerical
         * fluxes and sources (e.g., they have been multiplied by dt) so the hyperbolic flux
         * synchronization routine should be consistent with that.
         */
        virtual void
        synchronizeFluxes(
            hier::Patch& patch,
            const double time,
            const double dt) = 0;
        
        /**
         * This is an optional routine for user to process any application-specific patch strategy
         * data BEFORE patches are advanced on the given level. This routine is called after patch
         * boundary data is filled (i.e., ghosts) and before computeFluxesAndSourcesOnPatch(). The
         * arguments are:
         * level          -- level that will be advanced,
         * current_time   -- current integration time,
         * dt             -- current time increment,
         * first_step     -- boolean flag that is true if advance is first in time step sequence on
         *                   level (i.e., previous advance step was on another level, false otherwise,
         * last_step      -- boolean flag that is true if advance is last in time step sequence on
         *                   level (i.e., synchronization with coarser level will occur immediately
         *                   after this advance),
         * regrid_advance -- boolean flag that is true if the advance is during a regridding phase
         *                   (i.e., the advance is not used to integrate data on the hierarchy) in
         *                   which case the results of the advance will be discarded.
         *
         * Note that when this routine is called, the scratch data is filled on all patches (i.e.,
         * ghost cells) and that data is the same as the current level data on all patch interiors.
         * That is, both scratch and current data correspond to current_time.
         *
         * Note that this function is not pure virtual. It is given a dummy implementation here so
         * that users may ignore it when inheriting from this class.
         */
        virtual void
        preprocessAdvanceLevelState(
            const boost::shared_ptr<hier::PatchLevel>& level,
            double current_time,
            double dt,
            bool first_step,
            bool last_step,
            bool regrid_advance);
        
        /**
         * This is an optional routine for user to process any application-specific patch strategy
         * data AFTER patches are advanced on the given level. This routine is called after looping
         * over advanceSingleStepOnPatch() is completed and before computeStableDtOnPatch(). The
         * arguments are:
         * level          -- level that will be advanced,
         * current_time   -- current integration time,
         * dt             -- current time increment,
         * first_step     -- boolean flag that is true if advance is first in time step sequence on
         *                   level (i.e., previous advance step was on another level, false otherwise),
         * last_step      -- boolean flag that is true if advance is last in time step sequence on
         *                   level (i.e., synchronization with coarser level will occur immediately
         *                   after this advance),
         * regrid_advance -- boolean flag that is true if the advance is during a regridding phase
         *                   (i.e., the advance is not used to integrate data on the hierarchy) in
         *                   which case the results of the advance will be discarded.
         *
         * Note that when this routine is called, the scratch data is filled on all patches (i.e.,
         * ghost cells) and that data is the same as the new level data on all patch interiors. That
         * is, both scratch and new data correspond to current_time + dt on patch interiors. The
         * current data and ghost values correspond to the current_time.
         *
         * Note that this function is not pure virtual. It is given a dummy implementation here so
         * that users may ignore it when inheriting from this class.
         */
        virtual void
        postprocessAdvanceLevelState(
            const boost::shared_ptr<hier::PatchLevel>& level,
            double current_time,
            double dt,
            bool first_step,
            bool last_step,
            bool regrid_advance);
        
        /**
         * This is an optional routine for user to process any application-specific patch strategy
         * data BEFORE cells are tagged on the given level using value detector.
         */
        virtual void
        preprocessTagCellsValueDetector(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const double regrid_time,
            const bool initial_error,
            const bool uses_gradient_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * This is an optional routine for user to process any application-specific patch strategy
         * data AFTER cells are tagged on the given level using value detector.
         */
        virtual void
        postprocessTagCellsValueDetector(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const double regrid_time,
            const bool initial_error,
            const bool uses_gradient_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Tag cells on the given patch that require refinement based on application-specific numerical
         * quantities. The tag index argument indicates the index of the tag data on the patch data
         * array. The boolean argument initial_error is true if tagging is being done at the initial
         * simulation time; otherwise, it is false.
         *
         * The boolean uses_gradient_detector is true when gradient detector is used in additional
         * to the value detector. The boolearn uses_multiresolution_detector_too is true when
         * multiresolution detector is used in addition to the value detector. The boolean
         * uses_integral_detector_too is true when integral detector is used in addition to the value
         * detector. The boolean uses_richardson_extrapolation_too is true when Richardson
         * extrapolation is used in addition to the value detector. These flags help users manage
         * multiple regridding criteria.
         *
         * Note that this function is not pure virtual. It is given a dummy implementation here so
         * that users may ignore it when inheriting from this class.
         */
        virtual void
        tagCellsOnPatchValueDetector(
           hier::Patch& patch,
           const double regrid_time,
           const bool initial_error,
           const int tag_index,
           const bool uses_gradient_detector_too,
           const bool uses_multiresolution_detector_too,
           const bool uses_integral_detector_too,
           const bool uses_richardson_extrapolation_too);
        
        /**
         * This is an optional routine for user to process any application-specific patch strategy
         * data BEFORE cells are tagged on the given level using gradient detector.
         */
        virtual void
        preprocessTagCellsGradientDetector(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const double regrid_time,
            const bool initial_error,
            const bool uses_value_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * This is an optional routine for user to process any application-specific patch strategy
         * data AFTER cells are tagged on the given level using gradient detector.
         */
        virtual void
        postprocessTagCellsGradientDetector(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const double regrid_time,
            const bool initial_error,
            const bool uses_value_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Tag cells on the given patch that require refinement based on application-specific
         * numerical quantities. The tag index argument indicates the index of the tag data on the
         * patch data array. The boolean argument initial_error is true if tagging is being done at
         * the initial simulation time; otherwise, it is false.
         *
         * The boolean uses_value_detector is true when value detector is used in additional to the
         * gradient detector. The boolearn uses_multiresolution_detector_too is true when multiresolution
         * detector is used in addition to the gradient detector. The boolean uses_integral_detector_too
         * is true when integral detector is used in addition to the gradient detector. The boolean
         * uses_richardson_extrapolation_too is true when Richardson extrapolation is used in addition
         * to the gradient detector. These flags help users manage multiple regridding criteria.
         *
         * Note that this function is not pure virtual. It is given a dummy implementation here so
         * that users may ignore it when inheriting from this class.
         */
        virtual void
        tagCellsOnPatchGradientDetector(
           hier::Patch& patch,
           const double regrid_time,
           const bool initial_error,
           const int tag_index,
           const bool uses_value_detector_too,
           const bool uses_multiresolution_detector_too,
           const bool uses_integral_detector_too,
           const bool uses_richardson_extrapolation_too);
        
        /**
         * This is an optional routine for user to process any application-specific patch strategy
         * data BEFORE cells are tagged on the given level using multiresolution detector.
         */
        virtual void
        preprocessTagCellsMultiresolutionDetector(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const double regrid_time,
            const bool initial_error,
            const bool uses_value_detector_too,
            const bool uses_gradient_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * This is an optional routine for user to process any application-specific patch strategy
         * data AFTER cells are tagged on the given level using multiresolution detector.
         */
        virtual void
        postprocessTagCellsMultiresolutionDetector(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const double regrid_time,
            const bool initial_error,
            const bool uses_value_detector_too,
            const bool uses_gradient_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Tag cells on the given patch that require refinement based on application-specific numerical
         * quantities. The tag index argument indicates the index of the tag data on the patch data
         * array. The boolean argument initial_error is true if tagging is being done at the initial
         * simulation time; otherwise, it is false.
         *
         * The boolean uses_value_detector_too is true when value detector is used in addition to
         * the multiresolution detector. The boolean uses_gradient_detector_too is true when gradient
         * detector is used in addition to the multiresolution detector. The boolean
         * uses_integral_detector_too is true when integral detector is used in addition to the
         * multiresolution detector. The boolean uses_richardson_extrapolation_too is true when
         * Richardson extrapolation is used in addition to the gradient detector. These flags help
         * users manage multiple regridding criteria.
         *
         * Note that this function is not pure virtual. It is given a dummy implementation here so
         * that users may ignore it when inheriting from this class.
         */
        virtual void
        tagCellsOnPatchMultiresolutionDetector(
           hier::Patch& patch,
           const double regrid_time,
           const bool initial_error,
           const int tag_index,
           const bool uses_value_detector_too,
           const bool uses_gradient_detector_too,
           const bool uses_integral_detector_too,
           const bool uses_richardson_extrapolation_too);
        
        /**
         * This is an optional routine for user to process any application-specific patch strategy
         * data BEFORE cells are tagged on the given level using integral detector.
         */
        virtual void
        preprocessTagCellsIntegralDetector(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const double regrid_time,
            const bool initial_error,
            const bool uses_value_detector_too,
            const bool uses_gradient_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * This is an optional routine for user to process any application-specific patch strategy
         * data AFTER cells are tagged on the given level using integral detector.
         */
        virtual void
        postprocessTagCellsIntegralDetector(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const double regrid_time,
            const bool initial_error,
            const bool uses_value_detector_too,
            const bool uses_gradient_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Tag cells on the given patch that require refinement based on application-specific numerical
         * quantities. The tag index argument indicates the index of the tag data on the patch data
         * array. The boolean argument initial_error is true if tagging is being done at the initial
         * simulation time; otherwise, it is false.
         *
         * The boolean uses_value_detector_too is true when value detector is used in addition to
         * the integral detector. The boolean uses_gradient_detector_too is true when gradient
         * detector is used in addition to the integral detector. The boolean
         * uses_multiresolution_detector_too is true when multiresolution detector is used in addition
         * to the integral detector. The boolean uses_richardson_extrapolation_too is true when
         * Richardson extrapolation is used in addition to the integral detector. These flags help
         * users manage multiple regridding criteria.
         *
         * Note that this function is not pure virtual. It is given a dummy implementation here so
         * that users may ignore it when inheriting from this class.
         */
        virtual void
        tagCellsOnPatchIntegralDetector(
           hier::Patch& patch,
           const double regrid_time,
           const bool initial_error,
           const int tag_index,
           const bool uses_value_detector_too,
           const bool uses_gradient_detector_too,
           const bool uses_multiresolution_detector_too,
           const bool uses_richardson_extrapolation_too);
        
        /**
         * Tag cells based from differences computed in the Richardson extrapolation. The Richardson
         * extrapolation algorithm creates a coarsened version of some hierarchy patch level and
         * advances data in time on both the coarsened patch level and the hierarchy level. This
         * routine takes the data resulting from the advance on both the coarse and fine levels,
         * compares them, and tags cells according to the difference.
         * \verbatim
         *              (2)
         *      n+1 ^------->x finish      (1) advanced_coarse
         *          |        ^             (2) coarsened_fine
         *  time  n -        |
         *          ^        |(1)
         *          |        |
         *          <--------o start
         *        fine     coarse
         * \endverbatim
         *
         * The patch supplied to this routine is on the coarsened level. However, the error_level_number
         * corresponds to the actual hierarchy level from which it was coarsened. Data resides on
         * this patch in two contexts - ``advanced_coarse'' and ``coarsened_fine''. Advanced coarse
         * is data advanced on the coarsened version of the level, while coarsened fine is the data
         * advanced on the fine level and then coarsened to the coarse level. The regrid time and
         * the time increment are given for the actual hierarchy level. The error coarsen ratio
         * argument is the ratio between the index spaces on the hierarchy level and the coarsened
         * hierarchy level. The boolean flag ``initial_error'' is true when the error estimation is
         * performed at the initial simulation time; i.e., when the hierarchy levels are being
         * constructed for the first time. The tag index argument is the index of the tag data on
         * the patch data array. The other boolean flag uses_gradient_detector_too is true when a
         * gradient detector scheme is used in addition to Richardson extrapolation. This flag helps
         * users manage multiple regridding criteria.
         *
         * The boolean uses_value_detector_too is true when value detector is used addition to the
         * Richardson extrapolation. The boolean uses_gradient_detector_too is true when gradient
         * detector is used in addition to the Richardson extrapolation. The boolean
         * uses_multiresolution_detector_too is true when multiresolution detector is used in addition
         * to the Richardson extrapolation. The boolean uses_integral_detector_too is true when
         * integral detector is used in addition to the Richardson extrapolation. These flags help
         * users manage multiple regridding criteria.
         * 
         * Note that this function is not pure virtual. It is given a dummy implementation here so
         * that users may ignore it when inheriting from this class.
         */
        virtual void
        tagCellsOnPatchRichardsonExtrapolation(
            hier::Patch& patch,
            const int error_level_number,
            const boost::shared_ptr<hier::VariableContext>& coarsened_fine,
            const boost::shared_ptr<hier::VariableContext>& advanced_coarse,
            const double regrid_time,
            const double deltat,
            const int error_coarsen_ratio,
            const bool initial_error,
            const int tag_index,
            const bool uses_value_detector_too,
            const bool uses_gradient_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_integral_detector_too);
        
        /**
         * Set user-defined boundary conditions at the physical domain boundary.
         */
        virtual void
        setPhysicalBoundaryConditions(
            hier::Patch& patch,
            const double fill_time,
            const hier::IntVector& ghost_width_to_fill) = 0;
        
        /**
         * Compute variables for computing the statistics of data.
         */
        virtual void
        computeStatisticsVariables(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy);
        
        /**
         * Filter variables for computing the statistics of data.
         */
        virtual void
        filterStatisticsVariables(
            const int level,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy);
        
        /**
         * Output the statistics of data.
         */
        virtual void
        outputDataStatistics(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time);
        
        /**
         * Return pointer to patch data context.
         */
        boost::shared_ptr<hier::VariableContext>
        getDataContext() const
        {
           return d_data_context;
        }
        
        /**
         * The Runge-Kutta integrator controls the context for the data to be used in the numerical
         * routines implemented in the concrete patch strategy. The setDataContext() allows the
         * integrator to set the context for data on a patch on which to operate.
         */
        void
        setDataContext(
           const boost::shared_ptr<hier::VariableContext>& context)
        {
           d_data_context = context;
        }
        
        /**
         * The clearDataContext() routine resets the data context to be null.
         */
        void
        clearDataContext()
        {
           d_data_context.reset();
        }
        
    private:
        boost::shared_ptr<hier::VariableContext> d_data_context;
        
};

#endif /* RUNGE_KUTTA_PATCH_STRATEGY_HPP */
