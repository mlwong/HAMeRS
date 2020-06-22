/*************************************************************************
 *
 * This file is modified from HyperbolicLevelIntegrator.c of the SAMRAI
 * distribution. For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Runge-Kutta integration routines for single level in AMR
 *                hierarchy
 *                (basic hyperbolic systems)
 *
 ************************************************************************/

#include "algs/integrator/RungeKuttaLevelIntegrator.hpp"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceDataFactory.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/OuterfaceData.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideDataFactory.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <cstdlib>
#include <fstream>
#include <string>

/*
 **************************************************************************************************
 *
 * External declarations for FORTRAN 77 routines used in flux synchronization process between
 * hierarchy levels.
 *
 **************************************************************************************************
 */
extern "C"
{
    
#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif
    
    // in upfluxsum1d.m4:
    void SAMRAI_F77_FUNC(upfluxsum1d, UPFLUXSUM1D) (
        const int&,
        const int&,
        const int&,
        const int&,
        const double *,
        double *);
    // in upfluxsum2d.m4:
    void SAMRAI_F77_FUNC(upfluxsumface2d0, UPFLUXSUMFACE2D0) (
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double *, double *);
    void SAMRAI_F77_FUNC(upfluxsumface2d1, UPFLUXSUMFACE2D1) (
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double *,
        double *);
    void SAMRAI_F77_FUNC(upfluxsumside2d0, UPFLUXSUMSIDE2D0) (
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double *,
        double *);
    void SAMRAI_F77_FUNC(upfluxsumside2d1, UPFLUXSUMSIDE2D1) (
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double *,
        double *);
    // in upfluxsum3d.m4:
    void SAMRAI_F77_FUNC(upfluxsumface3d0, UPFLUXSUMFACE3D0) (
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double *,
        double *);
    void SAMRAI_F77_FUNC(upfluxsumface3d1, UPFLUXSUMFACE3D1) (
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double *,
        double *);
    void SAMRAI_F77_FUNC(upfluxsumface3d2, UPFLUXSUMFACE3D2) (
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double *,
        double *);
    void SAMRAI_F77_FUNC(upfluxsumside3d0, UPFLUXSUMSIDE3D0) (
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double *,
        double *);
    void SAMRAI_F77_FUNC(upfluxsumside3d1, UPFLUXSUMSIDE3D1) (
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double *,
        double *);
    void SAMRAI_F77_FUNC(upfluxsumside3d2, UPFLUXSUMSIDE3D2) (
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const int&,
        const double *,
        double *);
}


const int RungeKuttaLevelIntegrator::RUNGE_KUTTA_LEVEL_INTEGRATOR_VERSION = 1;

bool RungeKuttaLevelIntegrator::s_barrier_after_error_bdry_fill_comm = true;

tbox::StartupShutdownManager::Handler
RungeKuttaLevelIntegrator::s_initialize_handler(
    RungeKuttaLevelIntegrator::initializeCallback,
    0,
    0,
    RungeKuttaLevelIntegrator::finalizeCallback,
    tbox::StartupShutdownManager::priorityTimers);

/*
 * Timers interspersed throughout the class.
 */
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_advance_bdry_fill_comm;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_error_bdry_fill_create;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_error_bdry_fill_comm;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_advance_mpi_reductions;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_get_level_dt_mpi_reductions;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_initialize_level_data;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_init_level_create_sched;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_init_level_fill_data;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_init_level_fill_interior;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_advance_bdry_fill_create;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_new_advance_bdry_fill_create;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_apply_value_detector;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_apply_gradient_detector;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_apply_multiresolution_detector;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_apply_integral_detector;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_tag_cells;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_coarsen_rich_extrap;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_get_level_dt;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_get_level_dt_sync;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_advance_level;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_new_advance_bdry_fill_comm;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_patch_num_kernel;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_advance_level_sync;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_std_level_sync;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_sync_new_levels;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_barrier_after_error_bdry_fill_comm;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_sync_initial_comm;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_sync_initial_create;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_coarsen_fluxsum_create;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_coarsen_fluxsum_comm;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_coarsen_sync_create;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_coarsen_sync_comm;
boost::shared_ptr<tbox::Timer> RungeKuttaLevelIntegrator::t_output_data_statistics;

#ifdef HLI_RECORD_STATS
/*
 * Statistics on number of cells and patches generated.
 */
std::vector<boost::shared_ptr<tbox::Statistic> > RungeKuttaLevelIntegrator::s_boxes_stat;
std::vector<boost::shared_ptr<tbox::Statistic> > RungeKuttaLevelIntegrator::s_cells_stat;
std::vector<boost::shared_ptr<tbox::Statistic> > RungeKuttaLevelIntegrator::s_timestamp_stat;
#endif

/*
 **************************************************************************************************
 *
 * This constructor sets the RungeKuttaPatchStrategy pointer and initializes integration parameters
 * to default values. Communication algorithms are created here too. Other data members are read
 * in from the input database or from the restart database corresponding to the specified object_name.
 *
 **************************************************************************************************
 */
RungeKuttaLevelIntegrator::RungeKuttaLevelIntegrator(
    const std::string& object_name,
    const boost::shared_ptr<tbox::Database>& input_db,
    RungeKuttaPatchStrategy* patch_strategy,
    bool use_time_refinement):
    d_patch_strategy(patch_strategy),
    d_object_name(object_name),
    d_use_time_refinement(use_time_refinement),
    d_use_cfl(true),
    d_cfl(tbox::MathUtilities<double>::getSignalingNaN()),
    d_cfl_init(tbox::MathUtilities<double>::getSignalingNaN()),
    d_lag_dt_computation(true),
    d_use_ghosts_for_dt(false),
    d_dt(tbox::MathUtilities<double>::getSignalingNaN()),
    d_flux_is_face(true),
    d_flux_face_registered(false),
    d_flux_side_registered(false),
    d_number_time_data_levels(2),
    d_scratch(hier::VariableDatabase::getDatabase()->getContext("SCRATCH")),
    d_current(hier::VariableDatabase::getDatabase()->getContext("CURRENT")),
    d_new(hier::VariableDatabase::getDatabase()->getContext("NEW")),
    d_plot_context(d_current),
    d_stats_tmp(hier::VariableDatabase::getDatabase()->getContext("STATS_TMP")),
    d_have_flux_on_level_zero(false),
    d_distinguish_mpi_reduction_costs(false)
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(patch_strategy != 0);
    
    tbox::RestartManager::getManager()->registerRestartItem(d_object_name,
        this);
    
    /*
     * Initialize object with data read from the input and restart databases.
     */
    bool from_restart = tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    getFromInput(input_db, from_restart);
    
    d_intermediate.resize(d_number_steps);
    d_intermediate_time_dep_data.resize(d_number_steps);
    d_intermediate_flux_var_data.resize(d_number_steps);
    d_intermediate_source_var_data.resize(d_number_steps);
    for (int sn = 0; sn < d_number_steps; sn++)
    {
        std::string context_string = "INTERMEDIATE_" + tbox::Utilities::intToString(sn);
        d_intermediate[sn] = hier::VariableDatabase::getDatabase()->getContext(context_string);
    }
    
    /*
     * If d_use_cfl is false, d_use_time_refinement is turned off automatically.
     */
    if (!d_use_cfl)
    {
        d_use_time_refinement = false;
    }
}


/*
 **************************************************************************************************
 *
 * Destructor tells the tbox::RestartManager to remove this object from the list of restart items.
 *
 **************************************************************************************************
 */
RungeKuttaLevelIntegrator::~RungeKuttaLevelIntegrator()
{
    tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
}


/*
 **************************************************************************************************
 *
 * Initialize integration data on all patches on level. This process is used at the start of the
 * simulation to set the initial hierarchy data and after adaptive regridding. In the second case,
 * the old level pointer points to the level that existed in the hierarchy before regridding. This
 * pointer may be null, in which case it is ignored. If it is non-null, then data is copied from the
 * old level to the new level before the old level is discarded.
 *
 * Note that we also allocate flux storage for the coarsest AMR hierarchy level here (i.e., level 0).
 * The time step sequence on level 0 is dictated by the user code; so to avoid any memory management
 * errors, flux storage on level 0 persists as long as the level does.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::initializeLevelData(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const boost::shared_ptr<hier::PatchLevel>& old_level,
    const bool allocate_data)
{
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(!old_level || level_number == old_level->getLevelNumber());
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    
    boost::shared_ptr<hier::PatchLevel> level(
        hierarchy->getPatchLevel(level_number));
    
    const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());
    mpi.Barrier();
    t_initialize_level_data->start();
    
    /*
     * Allocate storage needed to initialize level and fill data
     * from coarser levels in AMR hierarchy, potentially. Since
     * time gets set when we allocate data, re-stamp it to current
     * time if we don't need to allocate.
     */
    if (allocate_data)
    {
        level->allocatePatchData(d_new_patch_init_data, init_data_time);
        level->allocatePatchData(d_old_time_dep_data, init_data_time);
    }
    else
    {
        level->setTime(init_data_time, d_new_patch_init_data);
    }
    
    /*
     * Create schedules for filling new level and fill data.
     */
    level->getBoxLevel()->getMPI().Barrier();
    
    if ((level_number > 0) || old_level)
    {
        t_init_level_create_sched->start();
        
        boost::shared_ptr<xfer::RefineSchedule> sched(
            d_fill_new_level->createSchedule(
                level,
                old_level,
                level_number - 1,
                hierarchy,
                d_patch_strategy));
        
        mpi.Barrier();
        
        t_init_level_create_sched->stop();
        
        d_patch_strategy->setDataContext(d_scratch);
        
        t_init_level_fill_data->start();
        sched->fillData(init_data_time);
        mpi.Barrier();
        t_init_level_fill_data->stop();
        
        d_patch_strategy->clearDataContext();
    }
    
    if ((d_number_time_data_levels == 3) && can_be_refined)
    {
        hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
     
        for (hier::PatchLevel::iterator ip(level->begin());
             ip != level->end();
             ip++)
        {
            const boost::shared_ptr<hier::Patch>& patch = *ip;
            
            std::list<boost::shared_ptr<hier::Variable> >::iterator time_dep_var =
               d_time_dep_variables.begin();
            
            while (time_dep_var != d_time_dep_variables.end())
            {
                int old_indx = variable_db->mapVariableAndContextToIndex(
                    *time_dep_var,
                    d_old);
                int cur_indx = variable_db->mapVariableAndContextToIndex(
                    *time_dep_var,
                    d_current);
                
                patch->setPatchData(old_indx, patch->getPatchData(cur_indx));
                
                time_dep_var++;
            }
        }
    }

    mpi.Barrier();
    t_init_level_fill_interior->start();
    
    /*
     * Initialize data on patch interiors.
     */
    d_patch_strategy->setDataContext(d_current);
    
    for (hier::PatchLevel::iterator p(level->begin()); p != level->end(); p++)
    {
        const boost::shared_ptr<hier::Patch>& patch = *p;
        
        patch->allocatePatchData(d_temp_var_scratch_data, init_data_time);
        
        d_patch_strategy->initializeDataOnPatch(
            *patch,
            init_data_time,
            initial_time);
        
        patch->deallocatePatchData(d_temp_var_scratch_data);
    }
    d_patch_strategy->clearDataContext();
    mpi.Barrier();
    t_init_level_fill_interior->stop();
    
    t_initialize_level_data->stop();
}


/*
 **************************************************************************************************
 *
 * Reset hierarchy configuration information where the range of new hierarchy levels is specified.
 * The information updated involves the cached communication schedules maintained by the algorithm.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::resetHierarchyConfiguration(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    NULL_USE(finest_level);
    
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((coarsest_level >= 0)
        && (coarsest_level <= finest_level)
        && (finest_level <= hierarchy->getFinestLevelNumber()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    for (int ln0 = 0; ln0 <= finest_level; ln0++)
    {
       TBOX_ASSERT(hierarchy->getPatchLevel(ln0));
    }
#endif
    
    int finest_hiera_level = hierarchy->getFinestLevelNumber();
    
    d_bdry_sched_advance.resize(finest_hiera_level + 1);
    d_bdry_sched_advance_new.resize(finest_hiera_level + 1);
    
    for (int ln = coarsest_level; ln <= finest_hiera_level; ln++)
    {
        boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
        
        t_advance_bdry_fill_create->start();
        
        d_bdry_sched_advance[ln] =
            d_bdry_fill_advance->createSchedule(
                level,
                ln - 1,
                hierarchy,
                d_patch_strategy,
                true);
        
        t_advance_bdry_fill_create->stop();
     
        if (!d_lag_dt_computation && d_use_ghosts_for_dt)
        {
            t_new_advance_bdry_fill_create->start();
            
            d_bdry_sched_advance_new[ln] =
                d_bdry_fill_advance_new->createSchedule(
                    level,
                    ln - 1,
                    hierarchy,
                    d_patch_strategy,
                    true);
            
            t_new_advance_bdry_fill_create->stop();
        }
    }
}


/*
 **************************************************************************************************
 *
 * Call patch routines to tag cells near large gradients. These cells will be refined.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::applyValueDetector(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) &&
                (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    
    t_apply_value_detector->start();
    
    boost::shared_ptr<hier::PatchLevel> level(
        hierarchy->getPatchLevel(level_number));
    
    level->allocatePatchData(d_saved_var_scratch_data, error_data_time);
    level->allocatePatchData(d_temp_var_scratch_data, error_data_time);
    
    d_patch_strategy->setDataContext(d_scratch);
    
    const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());
    
    t_error_bdry_fill_comm->start();
    d_bdry_sched_advance[level_number]->fillData(error_data_time);
    
    if (s_barrier_after_error_bdry_fill_comm)
    {
        t_barrier_after_error_bdry_fill_comm->start();
        mpi.Barrier();
        t_barrier_after_error_bdry_fill_comm->stop();
    }
    t_error_bdry_fill_comm->stop();
    
    d_patch_strategy->preprocessTagCellsValueDetector(
        hierarchy,
        level_number,
        error_data_time,
        initial_time,
        uses_gradient_detector_too,
        uses_multiresolution_detector_too,
        uses_integral_detector_too,
        uses_richardson_extrapolation_too);
    
    t_tag_cells->start();
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch>& patch = *ip;
        d_patch_strategy->tagCellsOnPatchValueDetector(
            *patch,
            error_data_time,
            initial_time,
            tag_index,
            uses_gradient_detector_too,
            uses_multiresolution_detector_too,
            uses_integral_detector_too,
            uses_richardson_extrapolation_too);
    }
    t_tag_cells->stop();
    
    d_patch_strategy->postprocessTagCellsValueDetector(
        hierarchy,
        level_number,
        error_data_time,
        initial_time,
        uses_gradient_detector_too,
        uses_multiresolution_detector_too,
        uses_integral_detector_too,
        uses_richardson_extrapolation_too);
    
    d_patch_strategy->clearDataContext();
    
    copyTimeDependentData(level, d_scratch, d_current);
    
    level->deallocatePatchData(d_temp_var_scratch_data);
    level->deallocatePatchData(d_saved_var_scratch_data);
    
    t_apply_value_detector->stop();
}


void
RungeKuttaLevelIntegrator::applyGradientDetector(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_value_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) &&
                (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    
    t_apply_gradient_detector->start();
    
    boost::shared_ptr<hier::PatchLevel> level(
        hierarchy->getPatchLevel(level_number));
    
    level->allocatePatchData(d_saved_var_scratch_data, error_data_time);
    level->allocatePatchData(d_temp_var_scratch_data, error_data_time);
    
    d_patch_strategy->setDataContext(d_scratch);
    
    const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());
    
    t_error_bdry_fill_comm->start();
    d_bdry_sched_advance[level_number]->fillData(error_data_time);
    
    if (s_barrier_after_error_bdry_fill_comm)
    {
        t_barrier_after_error_bdry_fill_comm->start();
        mpi.Barrier();
        t_barrier_after_error_bdry_fill_comm->stop();
    }
    t_error_bdry_fill_comm->stop();
    
    d_patch_strategy->preprocessTagCellsGradientDetector(
        hierarchy,
        level_number,
        error_data_time,
        initial_time,
        uses_value_detector_too,
        uses_multiresolution_detector_too,
        uses_integral_detector_too,
        uses_richardson_extrapolation_too);
    
    t_tag_cells->start();
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch>& patch = *ip;
        d_patch_strategy->tagCellsOnPatchGradientDetector(
            *patch,
            error_data_time,
            initial_time,
            tag_index,
            uses_value_detector_too,
            uses_multiresolution_detector_too,
            uses_integral_detector_too,
            uses_richardson_extrapolation_too);
    }
    t_tag_cells->stop();
    
    d_patch_strategy->postprocessTagCellsGradientDetector(
        hierarchy,
        level_number,
        error_data_time,
        initial_time,
        uses_value_detector_too,
        uses_multiresolution_detector_too,
        uses_integral_detector_too,
        uses_richardson_extrapolation_too);
    
    d_patch_strategy->clearDataContext();
    
    copyTimeDependentData(level, d_scratch, d_current);
    
    level->deallocatePatchData(d_temp_var_scratch_data);
    level->deallocatePatchData(d_saved_var_scratch_data);
    
    t_apply_gradient_detector->stop();
}


/*
 **************************************************************************************************
 *
 * Call patch routines to tag cells with multiresolution detector. These cells will be refined.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::applyMultiresolutionDetector(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) &&
                (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    
    t_apply_multiresolution_detector->start();
    
    boost::shared_ptr<hier::PatchLevel> level(
        hierarchy->getPatchLevel(level_number));
    
    level->allocatePatchData(d_saved_var_scratch_data, error_data_time);
    level->allocatePatchData(d_temp_var_scratch_data, error_data_time);
    
    d_patch_strategy->setDataContext(d_scratch);
    
    const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());
    
    t_error_bdry_fill_comm->start();
    d_bdry_sched_advance[level_number]->fillData(error_data_time);
    
    if (s_barrier_after_error_bdry_fill_comm)
    {
        t_barrier_after_error_bdry_fill_comm->start();
        mpi.Barrier();
        t_barrier_after_error_bdry_fill_comm->stop();
    }
    t_error_bdry_fill_comm->stop();
    
    d_patch_strategy->preprocessTagCellsMultiresolutionDetector(
        hierarchy,
        level_number,
        error_data_time,
        initial_time,
        uses_value_detector_too,
        uses_gradient_detector_too,
        uses_integral_detector_too,
        uses_richardson_extrapolation_too);
    
    t_tag_cells->start();
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch>& patch = *ip;
        d_patch_strategy->tagCellsOnPatchMultiresolutionDetector(
            *patch,
            error_data_time,
            initial_time,
            tag_index,
            uses_value_detector_too,
            uses_gradient_detector_too,
            uses_integral_detector_too,
            uses_richardson_extrapolation_too);
    }
    t_tag_cells->stop();
    
    d_patch_strategy->postprocessTagCellsMultiresolutionDetector(
        hierarchy,
        level_number,
        error_data_time,
        initial_time,
        uses_value_detector_too,
        uses_gradient_detector_too,
        uses_integral_detector_too,
        uses_richardson_extrapolation_too);
    
    d_patch_strategy->clearDataContext();
    
    copyTimeDependentData(level, d_scratch, d_current);
    
    level->deallocatePatchData(d_temp_var_scratch_data);
    level->deallocatePatchData(d_saved_var_scratch_data);
    
    t_apply_multiresolution_detector->stop();
}


/*
 **************************************************************************************************
 *
 * Call patch routines to tag cells with integral detector. These cells will be refined.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::applyIntegralDetector(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) &&
                (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    
    t_apply_integral_detector->start();
    
    boost::shared_ptr<hier::PatchLevel> level(
        hierarchy->getPatchLevel(level_number));
    
    level->allocatePatchData(d_saved_var_scratch_data, error_data_time);
    level->allocatePatchData(d_temp_var_scratch_data, error_data_time);
    
    d_patch_strategy->setDataContext(d_scratch);
    
    const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());
    
    t_error_bdry_fill_comm->start();
    d_bdry_sched_advance[level_number]->fillData(error_data_time);
    
    if (s_barrier_after_error_bdry_fill_comm)
    {
        t_barrier_after_error_bdry_fill_comm->start();
        mpi.Barrier();
        t_barrier_after_error_bdry_fill_comm->stop();
    }
    t_error_bdry_fill_comm->stop();
    
    d_patch_strategy->preprocessTagCellsIntegralDetector(
        hierarchy,
        level_number,
        error_data_time,
        initial_time,
        uses_value_detector_too,
        uses_gradient_detector_too,
        uses_multiresolution_detector_too,
        uses_richardson_extrapolation_too);
    
    t_tag_cells->start();
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch>& patch = *ip;
        d_patch_strategy->tagCellsOnPatchIntegralDetector(
            *patch,
            error_data_time,
            initial_time,
            tag_index,
            uses_value_detector_too,
            uses_gradient_detector_too,
            uses_multiresolution_detector_too,
            uses_richardson_extrapolation_too);
    }
    t_tag_cells->stop();
    
    d_patch_strategy->postprocessTagCellsIntegralDetector(
        hierarchy,
        level_number,
        error_data_time,
        initial_time,
        uses_value_detector_too,
        uses_gradient_detector_too,
        uses_multiresolution_detector_too,
        uses_richardson_extrapolation_too);
    
    d_patch_strategy->clearDataContext();
    
    copyTimeDependentData(level, d_scratch, d_current);
    
    level->deallocatePatchData(d_temp_var_scratch_data);
    level->deallocatePatchData(d_saved_var_scratch_data);
    
    t_apply_integral_detector->stop();
}


/*
 **************************************************************************************************
 *
 * Call patch routines to tag cells for refinement using Richardson extrapolation. Richardson
 * extrapolation requires two copies of the solution to compare. The NEW context holds the solution
 * computed on the fine level and coarsened, whereas the CURRENT context holds the solution integrated
 * on the coarse level after coarsening the initial data from the fine level.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::applyRichardsonExtrapolation(
    const boost::shared_ptr<hier::PatchLevel>& level,
    const double error_data_time,
    const int tag_index,
    const double deltat,
    const int error_coarsen_ratio,
    const bool initial_time,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too)
{
    TBOX_ASSERT(level);
    
    /*
     * Compare solutions computed on level (stored in NEW context) and on
     * the coarser level (stored in CURR context) on the patches of the
     * coarser level. The patch strategy implements the compare operations
     * performed on each patch.
     */
    
    const int error_level_number =
        level->getNextCoarserHierarchyLevelNumber() + 1;
    
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch>& patch = *ip;
        
        d_patch_strategy->tagCellsOnPatchRichardsonExtrapolation(
            *patch,
            error_level_number,
            d_new,                                     //  finer context
            d_current,                                 //  coarser context
            error_data_time,
            deltat,
            error_coarsen_ratio,
            initial_time,
            tag_index,
            uses_value_detector_too,
            uses_gradient_detector_too,
            uses_multiresolution_detector_too,
            uses_integral_detector_too);
    }
}


/*
 **************************************************************************************************
 *
 * The Richardson extrapolation algorithm requires a coarsened version of the level on which error
 * estiamtion is performed. This routine is used to coarsen data from a level in the AMR hierarchy
 * to some coarsened version of it. Note that this routine will be called twice. The before_advance
 * boolean argument indicates whether data is set on the coarse level by coarsening the "old" time
 * level solution or by coarsening the "new" solution on the fine level (i.e., after it has been
 * advanced).
 *
 * The contexts used for coarsening old data depend on the number of time levels. We always want to
 * use data at the oldest time on the fine level, coarsened to the CURRENT context on the coarse
 * level. Thus, if the problem uses two time levels, we coarsen data from CURRENT on fine level
 * (since CURRENT is the oldest time maintained) to CURRENT on the coarse level. If the problem uses
 * three time levels, we coarsen from OLD on the fine level (since OLD is the time maintained) to
 * CURRENT on the coarse level.
 *
 * When the boolean is false, indicating we are operating at the new time, we coarsen the time
 * advanced solution at the NEW context on the fine level to the NEW context on the coarse level so
 * that they may be compared later.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::coarsenDataForRichardsonExtrapolation(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const boost::shared_ptr<hier::PatchLevel>& coarse_level,
    const double coarsen_data_time,
    const bool before_advance)
{
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
    TBOX_ASSERT(coarse_level);
    TBOX_ASSERT_OBJDIM_EQUALITY2(*hierarchy, *coarse_level);
    
    t_coarsen_rich_extrap->start();
    
    const hier::IntVector& zero_vector(hier::IntVector::getZero(hierarchy->getDim()));
    
    boost::shared_ptr<hier::PatchLevel> hier_level(
        hierarchy->getPatchLevel(level_number));
    
    hier::IntVector coarsen_ratio(hierarchy->getDim());
    if (coarse_level->getRatioToLevelZero() < zero_vector)
    {
        if (hier_level->getRatioToLevelZero() < zero_vector)
        {
            coarsen_ratio = coarse_level->getRatioToLevelZero()
                / hier_level->getRatioToLevelZero();
            TBOX_ASSERT(
                coarsen_ratio * hier_level->getRatioToLevelZero() ==
                coarse_level->getRatioToLevelZero());
        }
        else
        {
            coarsen_ratio = -coarse_level->getRatioToLevelZero()
                / hier_level->getRatioToLevelZero();
            TBOX_ASSERT(
                -coarsen_ratio * hier_level->getRatioToLevelZero() ==
                coarse_level->getRatioToLevelZero());
        }
    }
    else
    {
        coarsen_ratio = hier_level->getRatioToLevelZero()
            / coarse_level->getRatioToLevelZero();
        TBOX_ASSERT(
            coarsen_ratio * coarse_level->getRatioToLevelZero() ==
            hier_level->getRatioToLevelZero());
    }
    
    /*
     * Compute the width needed for Connectors. The peer width for coarse<==>fine can be equivalent
     * to the width for fine<==>fine in the hierarcy, because the coarse level is just the coarsened
     * fine level. We just have to convert the width to the correct refinement ratio before initializing
     * the Connectors.
     */
    const hier::IntVector peer_connector_width =
        hierarchy->getRequiredConnectorWidth(
            level_number,
            level_number,
            true);
    
    const hier::IntVector c_to_f_width =
       hier::IntVector::ceilingDivide(peer_connector_width, coarsen_ratio);
    
    const hier::IntVector f_to_c_width(c_to_f_width * coarsen_ratio);
    
    coarse_level->findConnectorWithTranspose(
        *hier_level,
        c_to_f_width,
        f_to_c_width,
        hier::CONNECTOR_CREATE);
    
    if (before_advance)
    {
        coarse_level->allocatePatchData(
            d_new_patch_init_data,
            coarsen_data_time);
        
        if (d_number_time_data_levels == 3)
        {
            d_patch_strategy->setDataContext(d_old);
        }
        else
        {
            d_patch_strategy->setDataContext(d_current);
        }
     
        d_coarsen_rich_extrap_init->createSchedule(
            coarse_level,
            hier_level,
            d_patch_strategy)->coarsenData();
        
        d_patch_strategy->clearDataContext();
    }
    else
    {
        coarse_level->allocatePatchData(
            d_new_time_dep_data,
            coarsen_data_time);
     
        d_patch_strategy->setDataContext(d_new);
     
        d_coarsen_rich_extrap_final->createSchedule(
            coarse_level,
            hier_level,
            d_patch_strategy)
                ->coarsenData();
        
        d_patch_strategy->clearDataContext();
    }
    
    t_coarsen_rich_extrap->stop();
}


/*
 **************************************************************************************************
 *
 * Initialize level integrator by:
 *
 *   (1) Setting the number of time data levels based on needs of the gridding algorithm.
 *   (2) Invoking variable registration in patch strategy.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::initializeLevelIntegrator(
    const boost::shared_ptr<mesh::GriddingAlgorithmStrategy>& gridding_alg_strategy)
{
    d_gridding_alg =
        BOOST_CAST<mesh::GriddingAlgorithm, mesh::GriddingAlgorithmStrategy>(
            gridding_alg_strategy);
    
    TBOX_ASSERT(d_gridding_alg);
    
    d_number_time_data_levels = 2;
    
    if ((d_gridding_alg->getTagAndInitializeStrategy()->getErrorCoarsenRatio() < 1) ||
        (d_gridding_alg->getTagAndInitializeStrategy()->getErrorCoarsenRatio() > 3))
    {
        TBOX_ERROR("RungeKuttaLevelIntegrator::initializeLevelIntegrator "
                   << "error...\n"
                   << "   object name = "
                   << d_object_name
                   << "   gridding algorithm has bad error coarsen ratio"
                   << std::endl);
    }
    
    if ((d_gridding_alg->getTagAndInitializeStrategy()->everUsesTimeIntegration()) &&
        (d_gridding_alg->getTagAndInitializeStrategy()->getErrorCoarsenRatio() == 3))
    {
        d_number_time_data_levels = 3;
        d_old = hier::VariableDatabase::getDatabase()->getContext("OLD");
    }
    
    d_patch_strategy->registerModelVariables(this);
    
    d_patch_strategy->setupLoadBalancer(this, d_gridding_alg.get());
}


/*
 **************************************************************************************************
 *
 * Invoke dt calculation routines in patch strategy and take a min over all patches on the level.
 * The result will be the max of the next timestep on the level. If the boolean recompute_dt is true,
 * the max timestep on the level will be computed. If it is false, the method will simply access the
 * latest dt stored in the time refinement integrator.
 *
 **************************************************************************************************
 */
double
RungeKuttaLevelIntegrator::getLevelDt(
    const boost::shared_ptr<hier::PatchLevel>& level,
    const double dt_time,
    const bool initial_time)
{
    TBOX_ASSERT(level);
    
    double global_dt = 0.0;
    
    t_get_level_dt->start();
    
    if (d_use_cfl)
    {
        const tbox::SAMRAI_MPI& mpi(level->getBoxLevel()->getMPI());
        
        double dt = tbox::MathUtilities<double>::getMax();
        
        if (!d_use_ghosts_for_dt)
        {
            //tbox::plog << "!use ghosts for dt" << std::endl;
            
            d_patch_strategy->setDataContext(d_current);
            
            for (hier::PatchLevel::iterator p(level->begin());
                 p != level->end();
                 p++)
            {
                const boost::shared_ptr<hier::Patch>& patch = *p;
                
                patch->allocatePatchData(d_temp_var_scratch_data, dt_time);
                
                double patch_dt;
                patch_dt = d_patch_strategy->
                    computeStableDtOnPatch(
                        *patch,
                        initial_time,
                        dt_time);
                
                dt = tbox::MathUtilities<double>::Min(dt, patch_dt);
                //tbox::plog.precision(12);
                //tbox::plog << "Level " << level->getLevelNumber()
                //           << " Patch " << *p
                //           << " box " << patch->getBox()
                //           << " has patch_dt " << patch_dt
                //           << " dt " << dt
                //           << std::endl;
                
                patch->deallocatePatchData(d_temp_var_scratch_data);
            }
         
            d_patch_strategy->clearDataContext();
        }
        else
        {
            //tbox::plog << "use ghosts for dt" << std::endl;
            
            level->allocatePatchData(d_saved_var_scratch_data, dt_time);
            
            d_patch_strategy->setDataContext(d_scratch);
            
            t_advance_bdry_fill_comm->start();
            d_bdry_sched_advance[level->getLevelNumber()]->fillData(dt_time);
            t_advance_bdry_fill_comm->stop();
            
            for (hier::PatchLevel::iterator ip(level->begin());
                 ip != level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch>& patch = *ip;
                
                patch->allocatePatchData(d_temp_var_scratch_data, dt_time);
                
                double patch_dt;
                patch_dt = d_patch_strategy->
                    computeStableDtOnPatch(
                        *patch,
                        initial_time,
                        dt_time);
                
                dt = tbox::MathUtilities<double>::Min(dt, patch_dt);
                //tbox::plog.precision(12);
                //tbox::plog << "Level " << level->getLevelNumber()
                //           << " Patch " << *ip
                //           << " box " << patch->getBox()
                //           << " has patch_dt " << patch_dt
                //           << " dt " << dt
                //           << std::endl;
                
                patch->deallocatePatchData(d_temp_var_scratch_data);
            }
            
            d_patch_strategy->clearDataContext();
            
            /*
             * Copy data from scratch to current and de-allocate scratch storage. This may be excessive
             * here, but seems necessary if the computation of dt affects the state of the problem
             * solution. Also, this getLevelDt() routine is called at initialization only in most cases.
             */
         
            copyTimeDependentData(level, d_scratch, d_current);
         
            level->deallocatePatchData(d_saved_var_scratch_data);
        }
        
        t_get_level_dt_sync->start();
        
        if (d_distinguish_mpi_reduction_costs)
        {
            mpi.Barrier();
            t_get_level_dt_sync->stop();
            t_get_level_dt_mpi_reductions->start();
        }
        
        /*
         * The level time increment is a global min over all patches.
         */
        
        global_dt = dt;
        if (mpi.getSize() > 1)
        {
            mpi.AllReduce(&global_dt, 1, MPI_MIN);
        }
        global_dt *= tbox::MathUtilities<double>::Min(d_cfl_init, d_cfl);
        
        if (d_distinguish_mpi_reduction_costs)
        {
            t_get_level_dt_mpi_reductions->stop();
        }
        else
        {
            t_get_level_dt_sync->stop();
        }
    }
    else
    {
        global_dt = d_dt;
    }
    
    t_get_level_dt->stop();
    
    return global_dt;
}


/*
 **************************************************************************************************
 *
 * For the standard explicit integration algorithm for conservation laws, the fine time increment
 * is the coarse increment divided by the maximum mesh ratio (independent of level number).
 *
 **************************************************************************************************
 */
double
RungeKuttaLevelIntegrator::getMaxFinerLevelDt(
    const int finer_level_number,
    const double coarse_dt,
    const hier::IntVector& ratio)
{
    NULL_USE(finer_level_number);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    for (int id = 0; id < ratio.getDim().getValue(); id++)
    {
        TBOX_ASSERT(ratio(id) > 0);
    }
#endif
    return coarse_dt / double(ratio.max());
}


/*
 **************************************************************************************************
 *
 * Integrate data on all patches in patch level from current time to new time (new_time) using a
 * single time step. Before the advance can occur, proper ghost cell information is obtained for all
 * patches on the level. Then, local patches are advanced sequentially in the loop over patches. The
 * details of the routine are as follows:
 *
 *  0) Allocate storage for new time level data. Also, allocate necessary FLUX, flux integral and
 *     SOURCE storage if needed. (i.e., if regrid_advance is false, first_step is true, and coarser
 *     or finer level than current level exists in hierarchy.)
 *
 *  1) Scratch space is filled so that, for each patch, interior data and ghost cell bdry data
 *     correspond to specified time.
 *
 *  1a) Call user routines to pre-process advance data, if needed.
 *
 *  2) Advance in time by looping through Runge-Kutta sub-steps. In each sub-step,
 *
 *  2a) Copy time-dependent data from scratch space to corresponding intermediate space.
 *      
 *  2b) Compute next intermediate fluxes in corresponding intermediate space using data on patch +
 *      ghosts at given time.
 *
 *  2c) Advance one Runge-Kutta sub-step in scratch space and accumulate fluxs from intermediate
 *      fluxes using suitable weights for flux correction in flux synchronization later.
 *
 *  3) At the end of the loop, interior data is advanced to time = new_time.
 *
 *  3a) Call user routines to post-process advance data, if needed.
 *
 *  4) Compute next stable time increment for subsequent level advances:
 *
 *     4a) If (d_lag_dt_computation == true)
 *         {
 *            DO NOT RECOMPUTE characteristic data after advancing data on patch. Use characteristic
 *            data corresponding to current time level, computed prior to flux computation in dt
 *            calculation.
 *            If (d_use_ghosts_for_dt == true)
 *               - Compute dt using data on patch + ghosts at time.
 *            Else
 *               - Compute dt using data on patch interior ONLY.
 *         }
 *
 *     4b) Copy data from scratch space patch interior to new data storage for patch (i.e., at time =
 *         new_time).
 *
 *     4a) If (d_lag_dt_computation == false)
 *         {
 *            RECOMPUTE characteristic data after advancing data on patch. Use characteristic data
 *            corresponding to new time level in dt calculation.
 *            If (d_use_ghosts_for_dt == true)
 *               - Refill scratch space with new interior patch data and ghost cell bdry data
 *                  correspond to new time. (NOTE: This requires a new boundary schedule.)
 *               - Compute dt using data on patch+ghosts at new_time.
 *            Else
 *               - Compute dt using data on patch interior ONLY. (using patch interior data at
 *                 new_time)
 *         }
 *
 *  5) If (ln > 0), update flux integrals by adding patch bdry FLUXes to flux sums.
 *
 * Important Notes:
 *    1) In order to advance finer levels (if they exist), both old and new data for each patch on
 *       the level must be maintained.
 *    2) If the timestep is the first in the timestep loop on the level (indicated by first_step),
 *       then time interpolation is unnecessary to fill ghost cells from the next coarser level.
 *    3) The new dt is not calculated if regrid_advance is true. If this is the case, it is assumed
 *       that the results of the advance and the timestep calculation will be discarded (e.g., during
 *       regridding, or initialization). Also, allocation and post-processing of FLUX/flux integral
 *       data is not performed in this case.
 *
 **************************************************************************************************
 */
double
RungeKuttaLevelIntegrator::advanceLevel(
    const boost::shared_ptr<hier::PatchLevel>& level,
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const double current_time,
    const double new_time,
    const bool first_step,
    const bool last_step,
    const bool regrid_advance)
{
    TBOX_ASSERT(level);
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT(current_time <= new_time);
    TBOX_ASSERT_OBJDIM_EQUALITY2(*level, *hierarchy);
    
// HLI_RECORD_STATS is defined in RungeKuttaLevelIntegrator.h
#ifdef HLI_RECORD_STATS
    recordStatistics(*level, current_time);
#endif
    
    t_advance_level->start();
    
    const int level_number = level->getLevelNumber();
    const double dt = new_time - current_time;
    
    /*
     * (1) Allocate data needed for advancing level.
     * (2) Generate temporary communication schedule to fill ghost cells, if needed.
     * (3) Fill ghost cell data.
     * (4) Process flux and source storage before the advance.
     */
    
    level->allocatePatchData(d_new_time_dep_data, new_time);
    level->allocatePatchData(d_saved_var_scratch_data, current_time);
    
    for (int sn = 0; sn < d_number_steps; sn++)
    {
        level->allocatePatchData(d_intermediate_time_dep_data[sn], current_time);
        level->allocatePatchData(d_intermediate_flux_var_data[sn], current_time);
        level->allocatePatchData(d_intermediate_source_var_data[sn], current_time);
    }
    
    boost::shared_ptr<xfer::RefineSchedule> fill_schedule;
    
    if (!level->inHierarchy())
    {
        t_error_bdry_fill_create->start();
        
        const hier::OverlapConnectorAlgorithm oca;
        
        const int coarser_ln = level->getNextCoarserHierarchyLevelNumber();
        
        if (coarser_ln < 0)
        {
            // Don't use coarser level in boundary fill.
            if (d_number_time_data_levels == 3)
            {
                fill_schedule =
                    d_bdry_fill_advance_old->createSchedule(
                        level,
                        d_patch_strategy);
            }
            else
            {
                fill_schedule =
                    d_bdry_fill_advance->createSchedule(
                        level,
                        d_patch_strategy);
            }
        }
        else
        {
            // Use coarser level in boundary fill.
            if (d_number_time_data_levels == 3)
            {
                fill_schedule =
                    d_bdry_fill_advance_old->createSchedule(
                        level,
                        coarser_ln,
                        hierarchy,
                        d_patch_strategy,
                        true);
            }
            else
            {
                fill_schedule =
                    d_bdry_fill_advance->createSchedule(
                        level,
                        coarser_ln,
                        hierarchy,
                        d_patch_strategy,
                        true);
            }
        }
        t_error_bdry_fill_create->stop();
    }
    else
    {
        fill_schedule = d_bdry_sched_advance[level_number];
    }
    
    d_patch_strategy->setDataContext(d_scratch);
    
    if (regrid_advance)
    {
        t_error_bdry_fill_comm->start();
    }
    else
    {
        t_advance_bdry_fill_comm->start();
    }
    
    fill_schedule->fillData(current_time);
    
    if (regrid_advance)
    {
        t_error_bdry_fill_comm->stop();
    }
    else
    {
        t_advance_bdry_fill_comm->stop();
    }
    
    d_patch_strategy->clearDataContext();
    fill_schedule.reset();
    
    preprocessFluxAndSourceData(
        level,
        current_time,
        new_time,
        regrid_advance,
        first_step,
        last_step);
    
    /*
     * (5) Call user-routine to pre-process state data, if needed.
     * (6) Initialize all fluxes with zero values
     * (6) Advance solution on all level patches (scratch storage).
     *     In looping over Runge-Kutta sub-steps,
     *     (6a) Copy data from scatch data to the intermediate data.
     *          Dirchlet boundary conditions are applied at the coarse-fine boundaries of patches.
     *     (6b) Compute intermediate fluxes of current step.
     *     (6c) Advance one Runge-Kutta sub-step and accumulate the intermediate flux to the total
     *          flux during this whole Runge-Kutta step. Time-independent intermediate data of next
     *          Runge-Kutta step is stored in scratch context.
     * (7) Copy new solution to from scratch to new storage.
     * (8) Call user-routine to post-process state data, if needed.
     */
    
    t_patch_num_kernel->start();
    
    d_patch_strategy->preprocessAdvanceLevelState(
        level,
        current_time,
        dt,
        first_step,
        last_step,
        regrid_advance);
    
    t_patch_num_kernel->stop();
    
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch>& patch = *ip;
        
        patch->allocatePatchData(d_temp_var_scratch_data, current_time);
        
        // Fill all fluxes with zero values.
        
        std::list<boost::shared_ptr<hier::Variable> >::iterator flux_var =
            d_flux_variables.begin();
        
        while (flux_var != d_flux_variables.end())
        {            
            if (d_flux_is_face)
            {
                boost::shared_ptr<pdat::FaceData<double> > flux_data(
                    BOOST_CAST<pdat::FaceData<double>, hier::PatchData>(
                        patch->getPatchData(*flux_var, d_scratch)));
                
                TBOX_ASSERT(flux_data);
                flux_data->fillAll(0.0);
            }
            else
            {
                boost::shared_ptr<pdat::SideData<double> > flux_data(
                    BOOST_CAST<pdat::SideData<double>, hier::PatchData>(
                        patch->getPatchData(*flux_var, d_scratch)));
                
                TBOX_ASSERT(flux_data);
                flux_data->fillAll(0.0);
            }
            
            flux_var++;
        }
        
        // Fill all sources with zero values.
        
        std::list<boost::shared_ptr<hier::Variable> >::iterator source_var =
            d_source_variables.begin();
        
        while (source_var != d_source_variables.end())
        {
            boost::shared_ptr<pdat::CellData<double> > source_data(
                BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(*source_var, d_scratch)));
            
            source_data->fillAll(0.0);
            
            source_var++;
        }
    }
    
    const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());
    for (int sn = 0; sn < d_number_steps; sn++)
    {
        d_patch_strategy->setDataContext(d_intermediate[sn]);
        
        // Copy scratch data to intermediate data corresponding to current step.
        copyTimeDependentData(level, d_scratch, d_intermediate[sn]);
        
        /*
         * Fill the ghost cell data for current intemediate data factory.
         */
        
        boost::shared_ptr<xfer::RefineSchedule> fill_schedule_intermediate;
        
        if (sn > 0)
        {
            fill_schedule_intermediate = 
                d_bdry_fill_intermediate[sn]->createSchedule(
                    level,
                    d_patch_strategy);
            
            if (regrid_advance)
            {
                t_error_bdry_fill_comm->start();
            }
            else
            {
                t_advance_bdry_fill_comm->start();
            }
            
            fill_schedule_intermediate->fillData(current_time);
            
            if (regrid_advance)
            {
                t_error_bdry_fill_comm->stop();
            }
            else
            {
                t_advance_bdry_fill_comm->stop();
            }
        }
        
        d_patch_strategy->setDataContext(d_scratch);
            
        for (hier::PatchLevel::iterator ip(level->begin());
             ip != level->end();
             ip++)
        {
            const boost::shared_ptr<hier::Patch>& patch = *ip;
            
            t_patch_num_kernel->start();
            
            // Compute flux corresponding to this sub-step.
            d_patch_strategy->computeFluxesAndSourcesOnPatch(
                *patch,
                current_time,
                dt,
                sn,
                d_intermediate[sn]);
            
            // Advance a Runge-Kutta sub-step.
            d_patch_strategy->advanceSingleStepOnPatch(
                *patch,
                current_time,
                dt,
                d_alpha[sn],
                d_beta[sn],
                d_gamma[sn],
                d_intermediate);
            
            t_patch_num_kernel->stop();
        }
        
        fill_schedule_intermediate.reset();
    }
    
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch>& patch = *ip;
        
        patch->deallocatePatchData(d_temp_var_scratch_data);
    }
    
    d_patch_strategy->clearDataContext();
    
    level->setTime(new_time, d_saved_var_scratch_data);
    level->setTime(new_time, d_flux_var_data);
    
    copyTimeDependentData(level, d_scratch, d_new);
    
    t_patch_num_kernel->start();
    d_patch_strategy->postprocessAdvanceLevelState(
        level,
        current_time,
        dt,
        first_step,
        last_step,
        regrid_advance);
    t_patch_num_kernel->stop();
    
    /*
    * (9) If the level advance is for regridding, we compute the next timestep:
    *
    * (a) If the dt computation is lagged (i.e., we use pre-advance data to compute timestep), we
    *     reset scratch space on patch interiors if needed. Then, we set the strategy context to
    *     current or scratch depending on whether ghost values are used to compute dt.
    * (b) If the dt computation is not lagged (i.e., we use advanced data to compute timestep), we
    *     refill scratch space, including ghost data with new solution values if needed. Then, we
    *     set the strategy context to new or scratch depending on whether ghost values are used to
    *     compute dt.
    * (c) Then, we loop over patches and compute the dt on each patch.
    */
    
    double dt_next = tbox::MathUtilities<double>::getMax();
    
    if (!regrid_advance)
    {
        if (d_use_cfl)
        {
            if (d_lag_dt_computation)
            {
                if (d_use_ghosts_for_dt)
                {
                    d_patch_strategy->setDataContext(d_scratch);
                    copyTimeDependentData(level, d_current, d_scratch);
                }
                else
                {
                    d_patch_strategy->setDataContext(d_current);
                }
            }
            else
            {
                if (d_use_ghosts_for_dt)
                {
                    if (!d_bdry_sched_advance_new[level_number])
                    {
                        TBOX_ERROR(d_object_name
                            << ":  "
                            << "Attempt to fill new ghost data for timestep "
                            << "computation, but schedule not defined."
                            << std::endl);
                    }
                    d_patch_strategy->setDataContext(d_scratch);
                    t_new_advance_bdry_fill_comm->start();
                    d_bdry_sched_advance_new[level_number]->fillData(new_time);
                    t_new_advance_bdry_fill_comm->stop();
                }
                else
                {
                    d_patch_strategy->setDataContext(d_new);
                }
            }
          
            for (hier::PatchLevel::iterator ip(level->begin());
                 ip != level->end();
                 ip++)
            {
                const boost::shared_ptr<hier::Patch>& patch = *ip;
                
                patch->allocatePatchData(d_temp_var_scratch_data, new_time);
                
                // "false" argument indicates "initial_time" is false.
                t_patch_num_kernel->start();
                double patch_dt = d_patch_strategy->computeStableDtOnPatch(
                    *patch,
                    false,
                    new_time);
                t_patch_num_kernel->stop();
                
                dt_next = tbox::MathUtilities<double>::Min(dt_next, patch_dt);
                
                patch->deallocatePatchData(d_temp_var_scratch_data);
            }
            d_patch_strategy->clearDataContext();
        }
        else
        {
            dt_next = d_dt;
        }
      
    } // !regrid_advance
    
    level->deallocatePatchData(d_saved_var_scratch_data);
    
    for (int sn = 0; sn < d_number_steps; sn++)
    {
      level->deallocatePatchData(d_intermediate_time_dep_data[sn]);
      level->deallocatePatchData(d_intermediate_flux_var_data[sn]);
      level->deallocatePatchData(d_intermediate_source_var_data[sn]);
    }
    
    postprocessFluxAndSourceData(
        level,
        regrid_advance,
        first_step,
        last_step);
    
    t_advance_level->stop();
    
    t_advance_level_sync->start();
    
    if (d_distinguish_mpi_reduction_costs)
    {
        hierarchy->getMPI().Barrier();
        t_advance_level_sync->stop();
        t_advance_mpi_reductions->start();
    }
    
    double next_dt = dt_next;
    
    if (mpi.getSize() > 1)
    {
        mpi.AllReduce(&next_dt, 1, MPI_MIN);
    }
    
    if (d_use_cfl)
    {
        next_dt *= d_cfl;
    }
    
    if (d_distinguish_mpi_reduction_costs)
    {
        t_advance_mpi_reductions->stop();
    }
    else
    {
        t_advance_level_sync->stop();
    }
    
    return next_dt;
}


/*
 **************************************************************************************************
 * 
 * Synchronize data between patch levels according to the standard flux correction algorithm.
 * 
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::standardLevelSynchronization(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const double old_time)
{
    TBOX_ASSERT(hierarchy);
    
    std::vector<double> old_times(finest_level - coarsest_level + 1);
    for (int i = coarsest_level; i <= finest_level; i++)
    {
        old_times[i] = old_time;
    }
    
    standardLevelSynchronization(
        hierarchy,
        coarsest_level,
        finest_level,
        sync_time,
        old_times);
}


void
RungeKuttaLevelIntegrator::standardLevelSynchronization(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const std::vector<double>& old_times)
{
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((coarsest_level >= 0)
                && (coarsest_level < finest_level)
                && (finest_level <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(static_cast<int>(old_times.size()) >= finest_level);
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    for (int ln = coarsest_level; ln < finest_level; ln++)
    {
        TBOX_ASSERT(hierarchy->getPatchLevel(ln));
        TBOX_ASSERT(sync_time >= old_times[ln]);
    }
#endif
    TBOX_ASSERT(hierarchy->getPatchLevel(finest_level));
    
    t_std_level_sync->start();
    
    for (int fine_ln = finest_level; fine_ln > coarsest_level; --fine_ln)
    {
        const int coarse_ln = fine_ln - 1;
        
        boost::shared_ptr<hier::PatchLevel> fine_level(
           hierarchy->getPatchLevel(fine_ln));
        
        boost::shared_ptr<hier::PatchLevel> coarse_level(
           hierarchy->getPatchLevel(coarse_ln));
        
        synchronizeLevelWithCoarser(
            fine_level,
            coarse_level,
            sync_time,
            old_times[coarse_ln]);
        
        fine_level->deallocatePatchData(d_fluxsum_data);
        fine_level->deallocatePatchData(d_flux_var_data);
        fine_level->deallocatePatchData(d_source_var_data);
        
        if (coarse_ln > coarsest_level)
        {
            coarse_level->deallocatePatchData(d_flux_var_data);
            coarse_level->deallocatePatchData(d_source_var_data);
        }
        else
        {
            if (coarsest_level == 0)
            {
                coarse_level->deallocatePatchData(d_flux_var_data);
                d_have_flux_on_level_zero = false;
                coarse_level->deallocatePatchData(d_source_var_data);
            }
        }
    }
    
    t_std_level_sync->stop();
}


/*
 **************************************************************************************************
 *                                                                       
 * Coarsen current solution data from finest hierarchy level specified down through the coarsest
 * hierarchy level specified, if initial_time is true (i.e., hierarchy is being constructed at
 * initial simulation time). After data is coarsened, the user's initialization routine is called
 * to reset data (as needed by the application) before that solution is further coarsened to the
 * next coarser level in the hierarchy. If initial_time is false, then this routine does nothing
 * In that case, interpolation of data from coarser levels is sufficient to set data on new levels
 * in the hierarchy during regridding.
 * 
 * NOTE: The fact that this routine does nothing when called at any time later than when the AMR
 *       hierarchy is constructed initially may need to change at some point based on application
 *       needs.
 * 
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::synchronizeNewLevels(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const bool initial_time)
{
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((coarsest_level >= 0)
                && (coarsest_level < finest_level)
                && (finest_level <= hierarchy->getFinestLevelNumber()));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    for (int ln = coarsest_level; ln <= finest_level; ln++)
    {
        TBOX_ASSERT(hierarchy->getPatchLevel(ln));
    }
#endif
    
    t_sync_new_levels->start();
    
    if (initial_time)
    {
        d_patch_strategy->setDataContext(d_current);
        
        for (int fine_ln = finest_level; fine_ln > coarsest_level; --fine_ln)
        {
            const int coarse_ln = fine_ln - 1;
            
            boost::shared_ptr<hier::PatchLevel> fine_level(
                hierarchy->getPatchLevel(fine_ln));
            
            boost::shared_ptr<hier::PatchLevel> coarse_level(
                hierarchy->getPatchLevel(coarse_ln));
            
            t_sync_initial_create->start();
            boost::shared_ptr<xfer::CoarsenSchedule> sched(
                d_sync_initial_data->createSchedule(
                    coarse_level,
                    fine_level,
                    d_patch_strategy));
            t_sync_initial_create->stop();
            
            t_sync_initial_comm->start();
            sched->coarsenData();
            t_sync_initial_comm->stop();
            
            for (hier::PatchLevel::iterator p(coarse_level->begin());
                 p != coarse_level->end();
                 p++)
            {
                const boost::shared_ptr<hier::Patch>& patch = *p;
                
                patch->allocatePatchData(d_temp_var_scratch_data, sync_time);
                
                d_patch_strategy->initializeDataOnPatch(
                    *patch,
                    sync_time,
                    initial_time);
                
                patch->deallocatePatchData(d_temp_var_scratch_data);
            }
        }
        
        d_patch_strategy->clearDataContext();
    } // if (initial_time)
    
    t_sync_new_levels->stop();
}


/*
 **************************************************************************************************
 *
 * Synchronize data between coarse and fine patch levels according to the standard flux correction
 * algorithm. The steps of the algorithm are:
 *
 *    (1) Replace coarse time-space flux integrals at coarse-fine boundaries with time-space flux
 *    integrals computed on fine level.
 *    (2) Repeat conservative difference on coarse level with corrected fluxes.
 *    (3) Conservatively coarsen solution on interior of fine level to coarse level.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::synchronizeLevelWithCoarser(
    const boost::shared_ptr<hier::PatchLevel>& fine_level,
    const boost::shared_ptr<hier::PatchLevel>& coarse_level,
    const double sync_time,
    const double coarse_sim_time)
{
    TBOX_ASSERT(fine_level);
    TBOX_ASSERT(coarse_level);
    TBOX_ASSERT(coarse_level->getLevelNumber() ==
       (fine_level->getLevelNumber() - 1));
    TBOX_ASSERT_OBJDIM_EQUALITY2(*fine_level, *coarse_level);
    TBOX_ASSERT(sync_time > coarse_sim_time);
    
    /*
     * Coarsen flux integrals around fine patch boundaries to coarser level and replace coarse flux
     * information where appropriate. NULL patch model is passed in to avoid over complicating coarsen
     * process; i.e. patch model is not needed in coarsening of flux integrals.
     */
    
    t_coarsen_fluxsum_create->start();
    boost::shared_ptr<xfer::CoarsenSchedule> sched(
        d_coarsen_fluxsum->createSchedule(
            coarse_level,
            fine_level,
            0));
    t_coarsen_fluxsum_create->stop();
    
    t_coarsen_fluxsum_comm->start();
    sched->coarsenData();
    t_coarsen_fluxsum_comm->stop();
    
    /*
     * Repeat conservative difference on coarser level.
     */
    coarse_level->allocatePatchData(d_saved_var_scratch_data, coarse_sim_time);
    coarse_level->setTime(coarse_sim_time, d_flux_var_data);
    
    d_patch_strategy->setDataContext(d_scratch);
    
    t_advance_bdry_fill_comm->start();
    d_bdry_sched_advance[coarse_level->getLevelNumber()]->
        fillData(coarse_sim_time);
    t_advance_bdry_fill_comm->stop();
    
    const double reflux_dt = sync_time - coarse_sim_time;
    
    for (hier::PatchLevel::iterator ip(coarse_level->begin());
         ip != coarse_level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch>& patch = *ip;
        
        patch->allocatePatchData(d_temp_var_scratch_data, coarse_sim_time);
            
        d_patch_strategy->synchronizeFluxes(*patch,
           coarse_sim_time,
           reflux_dt);
        patch->deallocatePatchData(d_temp_var_scratch_data);
    }
    
    d_patch_strategy->clearDataContext();
    
    copyTimeDependentData(coarse_level, d_scratch, d_new);
    
    coarse_level->deallocatePatchData(d_saved_var_scratch_data);
    
    /*
     * Coarsen time-dependent data from fine patch interiors to coarse patches.
     */
    
    t_coarsen_sync_create->start();
    sched = d_coarsen_sync_data->createSchedule(
        coarse_level,
        fine_level,
        d_patch_strategy);
    t_coarsen_sync_create->stop();
    
    d_patch_strategy->setDataContext(d_new);
    
    t_coarsen_sync_comm->start();
    sched->coarsenData();
    t_coarsen_sync_comm->stop();
    
    d_patch_strategy->clearDataContext();
}


/*
 **************************************************************************************************
 *
 * Reset time-dependent data on patch level by replacing current data with new. The boolean argument
 * is used for odd refinement ratios (in particular 3 used in certain applications).
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::resetTimeDependentData(
   const boost::shared_ptr<hier::PatchLevel>& level,
   const double new_time,
   const bool can_be_refined)
{
    TBOX_ASSERT(level);
    
    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
    
    double cur_time = 0.;
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch>& patch = *ip;
        
        std::list<boost::shared_ptr<hier::Variable> >::iterator time_dep_var =
            d_time_dep_variables.begin();
        
        while (time_dep_var != d_time_dep_variables.end())
        {
            int cur_indx =
                variable_db->mapVariableAndContextToIndex(
                    *time_dep_var,
                    d_current);
            int new_indx =
                variable_db->mapVariableAndContextToIndex(
                    *time_dep_var,
                    d_new);
            
            cur_time = patch->getPatchData(cur_indx)->getTime();
            
            if (can_be_refined && d_number_time_data_levels == 3)
            {
                int old_indx = variable_db->mapVariableAndContextToIndex(
                    *time_dep_var,
                    d_old);
                
                patch->setPatchData(old_indx, patch->getPatchData(cur_indx));
                
                patch->setPatchData(cur_indx, patch->getPatchData(new_indx));
            }
            else
            {
                if (d_number_time_data_levels == 3)
                {
                    int old_indx = variable_db->mapVariableAndContextToIndex(
                        *time_dep_var,
                        d_old);
                    
                    patch->setPatchData(old_indx, patch->getPatchData(cur_indx));
                }
                
                patch->setPatchData(cur_indx, patch->getPatchData(new_indx));
            }
            
            patch->deallocatePatchData(new_indx);
            
            time_dep_var++;
        }
    }
    
    level->setTime(new_time, d_new_patch_init_data);
    
    if (d_number_time_data_levels == 3)
    {
        level->setTime(cur_time, d_old_time_dep_data);
    }
}


/*
 **************************************************************************************************
 *
 * Discard new data on level. This is used primarily to reset patch data after error estimation
 * (e.g., Richardson extrapolation.)
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::resetDataToPreadvanceState(
    const boost::shared_ptr<hier::PatchLevel>& level)
{
    TBOX_ASSERT(level);
    
    /*
     * De-allocate new context
     */
    level->deallocatePatchData(d_new_time_dep_data);
}


/*
 **************************************************************************************************
 *
 * Register given variable with algorithm according to specified algorithm role (i.e., RK_VAR_TYPE).
 * Assignment of descriptor indices to variable lists, component selectors, and communication
 * algorithms takes place here. The different cases are:
 *
 * TIME_DEP:
 *             The number of factories depends on the number of time levels of data that must be
 *             stored on patches to satisfy regridding reqs. Currently, there are two possibilities:
 *
 *             (1) If the coarsen ratios between levels are even, the error coarsening ratio will
 *                 be two and so only two time levels of data must be maintained on every level but
 *                 the finest as usual.
 *
 *             (2) If the coarsen ratios between levels are three, and time integration is used
 *                 during regridding (e.g., Richardson extrapolation), then three time levels of data
 *                 must be maintained on every level but the finest so that error estimation can be
 *                 executed properly.
 *
 *             In case (1), three factories are needed:
 *                          SCRATCH, CURRENT, NEW.
 *             In case (2), four factories are needed:
 *                          SCRATCH, OLD, CURRENT, NEW.
 *
 *             SCRATCH index is added to d_saved_var_scratch_data. CURRENT index is added to
 *             d_new_patch_init_data. NEW index is added to d_new_time_dep_data.
 *
 * INPUT:
 *             Only one time level of data is maintained and once values are set on patches, they
 *             do not change in time.
 *
 *             Two factories are needed: SCRATCH, CURRENT.
 *
 *             SCRATCH index is added to d_saved_var_scratch_data.
 *             CURRENT index is added to d_new_patch_init_data.
 *
 * NO_FILL:
 *             Only one time level of data is stored and no scratch space is used. Data may be set
 *             and manipulated at will in user routines. Data (including ghost values) is never
 *             touched outside of user routines.
 *
 *             Two factories are needed: CURRENT, SCRATCH.
 *
 *             CURRENT index is added to d_new_patch_init_data.
 *             SCRATCH index is needed only for temporary work space to fill new patch levels.
 *
 * FLUX:
 *             One factory is needed: SCRATCH.
 *
 *             SCRATCH index is added to d_flux_var_data.
 *
 *             Additionally, a variable for flux integral data is created for each FLUX variable. It
 *             has a single factory, SCRATCH, which is added to d_fluxsum_data.
 *
 * SOURCE:
 *             One factory is needed: SCRATCH.
 *
 *             SCRATCH index is added to d_source_var_data.
 *
 * TEMPORARY:
 *             One factory is needed: SCRATCH.
 *             SCRATCH index is added to d_temp_var_scratch_data.
 *
 * STATISTICS:
 *             Two factories is needed: SCRATCH, STATS_TMP.
 *             SCRATCH index is added to d_stats_var_scratch_data.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::registerVariable(
    const boost::shared_ptr<hier::Variable>& var,
    const hier::IntVector& ghosts,
    const hier::IntVector& ghosts_intermediate,
    const RK_VAR_TYPE& RK_v_type,
    const boost::shared_ptr<hier::BaseGridGeometry>& transfer_geom,
    const std::string& coarsen_name,
    const std::string& refine_name)
{
    const tbox::Dimension dim(ghosts.getDim());
    
    TBOX_ASSERT(var);
    TBOX_ASSERT(transfer_geom);
    TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, *var);
    
    TBOX_ASSERT(ghosts >= ghosts_intermediate);
    TBOX_ASSERT(ghosts.getDim() == ghosts_intermediate.getDim());
    
    if (!d_bdry_fill_advance)
    {
        /*
         * One-time set-up for communication algorithms.
         * We wait until this point to do this because we need a dimension.
         */
        d_bdry_fill_advance.reset(new xfer::RefineAlgorithm());
        d_bdry_fill_advance_new.reset(new xfer::RefineAlgorithm());
        d_bdry_fill_advance_old.reset(new xfer::RefineAlgorithm());
        
        d_bdry_fill_intermediate.resize(d_number_steps);
        for (int sn = 0; sn < d_number_steps; sn++)
        {
            d_bdry_fill_intermediate[sn].reset(new xfer::RefineAlgorithm());
        }
        
        d_fill_new_level.reset(new xfer::RefineAlgorithm());
        d_coarsen_fluxsum.reset(new xfer::CoarsenAlgorithm(dim));
        d_coarsen_sync_data.reset(new xfer::CoarsenAlgorithm(dim));
        d_sync_initial_data.reset(new xfer::CoarsenAlgorithm(dim));
        
        d_coarsen_rich_extrap_init.reset(new xfer::CoarsenAlgorithm(dim));
        d_coarsen_rich_extrap_final.reset(new xfer::CoarsenAlgorithm(dim));
    }
    
    if (!d_fill_statistics)
    {
        d_fill_statistics.reset(new xfer::RefineAlgorithm());
    }
    
    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
    
    const hier::IntVector& zero_ghosts(hier::IntVector::getZero(dim));
    
    d_all_variables.push_back(var);
    
    switch (RK_v_type)
    {
        case TIME_DEP:
        {
            d_time_dep_variables.push_back(var);
            
            int cur_id = variable_db->registerVariableAndContext(
                var,
                d_current,
                zero_ghosts);
            int new_id = variable_db->registerVariableAndContext(
                var,
                d_new,
                zero_ghosts);
            int scr_id = variable_db->registerVariableAndContext(
                var,
                d_scratch,
                ghosts);
            
            std::vector<int> intermediate_id(d_number_steps);
            for (int sn = 0; sn < d_number_steps; sn++)
            {
                intermediate_id[sn] = variable_db->registerVariableAndContext(
                    var,
                    d_intermediate[sn],
                    ghosts_intermediate);
            }
            
            d_saved_var_scratch_data.setFlag(scr_id);
            
            d_new_patch_init_data.setFlag(cur_id);
            
            d_new_time_dep_data.setFlag(new_id);
            
            for (int sn = 0; sn < d_number_steps; sn++)
            {
                d_intermediate_time_dep_data[sn].setFlag(intermediate_id[sn]);
            }
            
            /*
             * Register variable and context needed for restart.
             */
            hier::PatchDataRestartManager::getManager()->
                registerPatchDataForRestart(cur_id);
            
            /*
             * Set boundary fill schedules for time-dependent variable. If time interpolation operator
             * is non-NULL, regular advance bdry fill algorithm will time interpolate between current
             * and new data on coarser levels, and fill from current data on same level. New advance
             * bdry fill algorithm will time interpolate between current and new data on coarser
             * levels, and fill from new data on same level. If time interpolation operator is NULL,
             * regular and new bdry fill algorithms will use current and new data, respectively.
             */
            boost::shared_ptr<hier::RefineOperator> refine_op(
                transfer_geom->lookupRefineOperator(var, refine_name));
            
            boost::shared_ptr<hier::TimeInterpolateOperator> time_int(
                transfer_geom->lookupTimeInterpolateOperator(var));
            
            d_bdry_fill_advance->registerRefine(
                scr_id, cur_id, cur_id, new_id, scr_id, refine_op, time_int);
            d_bdry_fill_advance_new->registerRefine(
                scr_id, new_id, cur_id, new_id, scr_id, refine_op, time_int);
            d_fill_new_level->registerRefine(
                cur_id, cur_id, cur_id, new_id, scr_id, refine_op, time_int);
            
            /*
             * Set boundary fill schedules for data used in the intermediate steps of the Runge-Kutta
             * integration.
             */
            for (int sn = 0; sn < d_number_steps; sn++)
            {
                d_bdry_fill_intermediate[sn]->registerRefine(
                    intermediate_id[sn],
                    scr_id,
                    intermediate_id[sn],
                    refine_op);
            }
            
            /*
             * For data synchronization between levels, the coarsen algorithm will coarsen new data
             * on finer level to new data on coarser. Recall that coarser level data pointers will
             * not be reset until after synchronization so we always coarsen to new (see
             * synchronizeLevelWithCoarser routine).
             */
            
            boost::shared_ptr<hier::CoarsenOperator> coarsen_op(
                transfer_geom->lookupCoarsenOperator(var, coarsen_name));
            
            d_coarsen_sync_data->registerCoarsen(new_id, new_id, coarsen_op);
            
            d_sync_initial_data->registerCoarsen(cur_id, cur_id, coarsen_op);
            
            /*
             * Coarsen operations used in Richardson extrapolation. The init initializes data on
             * coarser level, before the coarse level advance. If two time levels are used, coarsening
             * occurs between the CURRENT context on both levels. If three levels are used, coarsening
             * occurs between the OLD context on the fine level and the CURRENT context on the coarse
             * level. The final coarsen algorithm coarsens data after it has been advanced on the
             * fine level to the NEW context on the coarser level.
             */
            
            if (d_number_time_data_levels == 3)
            {
                int old_id = variable_db->registerVariableAndContext(
                    var,
                    d_old,
                    zero_ghosts);
                d_old_time_dep_data.setFlag(old_id);
                
                d_bdry_fill_advance_old->registerRefine(
                    scr_id, cur_id, old_id, new_id, scr_id, refine_op, time_int);
                
                d_coarsen_rich_extrap_init->registerCoarsen(cur_id, old_id, coarsen_op);
            }
            else
            {
                d_coarsen_rich_extrap_init->registerCoarsen(cur_id, cur_id, coarsen_op);
            }
            
            d_coarsen_rich_extrap_final->registerCoarsen(new_id, new_id, coarsen_op);
            
            break;
        }
        case INPUT:
        {
            int cur_id = variable_db->registerVariableAndContext(
                var,
                d_current,
                zero_ghosts);
            int scr_id = variable_db->registerVariableAndContext(
                var,
                d_scratch,
                ghosts);
            
            d_saved_var_scratch_data.setFlag(scr_id);
            
            d_new_patch_init_data.setFlag(cur_id);
            
            /*
             * Register variable and context needed for restart.
             */
            hier::PatchDataRestartManager::getManager()->
            registerPatchDataForRestart(cur_id);
    
            /*
             * Bdry algorithms for input variables will fill from current only.
             */
            boost::shared_ptr<hier::RefineOperator> refine_op(
                transfer_geom->lookupRefineOperator(var, refine_name));
            
            d_bdry_fill_advance->registerRefine(
                scr_id, cur_id, scr_id, refine_op);
            d_bdry_fill_advance_new->registerRefine(
                scr_id, cur_id, scr_id, refine_op);
            d_fill_new_level->registerRefine(
                cur_id, cur_id, scr_id, refine_op);
            
            /*
             * At initialization, it may be necessary to coarsen INPUT data up through the hierarchy
             * so that all levels are consistent.
             */
            
            boost::shared_ptr<hier::CoarsenOperator> coarsen_op(
                transfer_geom->lookupCoarsenOperator(var, coarsen_name));
            
            d_sync_initial_data->registerCoarsen(cur_id, cur_id, coarsen_op);
            
            /*
             * Coarsen operation for setting initial data on coarser level in the Richardson
             * extrapolation algorithm.
             */
            
            d_coarsen_rich_extrap_init->
            registerCoarsen(cur_id, cur_id, coarsen_op);
            
            break;
        }
        case NO_FILL:
        {
            int cur_id = variable_db->registerVariableAndContext(
                var,
                d_current,
                ghosts);
            
            int scr_id = variable_db->registerVariableAndContext(
                var,
                d_scratch,
                ghosts);
            
            d_new_patch_init_data.setFlag(cur_id);
            
            /*
             * Register variable and context needed for restart.
             */
            hier::PatchDataRestartManager::getManager()->
                registerPatchDataForRestart(cur_id);
            
            boost::shared_ptr<hier::RefineOperator> refine_op(
                transfer_geom->lookupRefineOperator(var, refine_name));
            
            d_fill_new_level->registerRefine(
                cur_id, cur_id, scr_id, refine_op);
            
            /*
             * Coarsen operation for setting initial data on coarser level
             * in the Richardson extrapolation algorithm.
             */
            
            boost::shared_ptr<hier::CoarsenOperator> coarsen_op(
                transfer_geom->lookupCoarsenOperator(var, coarsen_name));
            
            d_coarsen_rich_extrap_init->registerCoarsen(cur_id, cur_id, coarsen_op);
            
            break;
        }
        case FLUX:
        {
            /*
             * Note that we force all flux variables to hold double precision data and be face- or
             * side-centered. Also, for each flux variable, a corresponding "fluxsum" variable is
             * created to manage synchronization of data betweeen patch levels in the hierarchy.
             */
            const boost::shared_ptr<pdat::FaceVariable<double> > face_var(
                boost::dynamic_pointer_cast<pdat::FaceVariable<double>,
                hier::Variable>(var));
            
            const boost::shared_ptr<pdat::SideVariable<double> > side_var(
                boost::dynamic_pointer_cast<pdat::SideVariable<double>,
                hier::Variable>(var));
            
            if (face_var)
            {
                if (d_flux_side_registered)
                {
                    TBOX_ERROR(d_object_name
                               << ":  "
                               << "Attempt to register FaceVariable when "
                               << "SideVariable already registered."
                               << std::endl);
                }
                
                d_flux_is_face = true;
            }
            else if(side_var)
            {
                if (d_flux_face_registered)
                {
                    TBOX_ERROR(d_object_name
                               << ":  "
                               << "Attempt to register SideVariable when "
                               << "FaceVariable already registered."
                               << std::endl);
                }
                
               d_flux_is_face = false;
            }
            else
            {
               TBOX_ERROR(d_object_name
                          << ":  "
                          << "Flux is neither face- or side-centered."
                          << std::endl);
            }
            
            d_flux_variables.push_back(var);
            
            int scr_id = variable_db->registerVariableAndContext(
                var,
                d_scratch,
                ghosts);
            
            d_flux_var_data.setFlag(scr_id);
            
            std::vector<int> intermediate_id(d_number_steps);
            for (int sn = 0; sn < d_number_steps; sn++)
            {
                intermediate_id[sn] = variable_db->registerVariableAndContext(
                    var,
                    d_intermediate[sn],
                    ghosts_intermediate);
            }
            
            for (int sn = 0; sn < d_number_steps; sn++)
            {
                d_intermediate_flux_var_data[sn].setFlag(intermediate_id[sn]);
            }
                
            std::string var_name = var->getName();
            std::string fs_suffix = "_fluxsum";
            std::string fsum_name = var_name;
            fsum_name += fs_suffix;
            
            boost::shared_ptr<hier::Variable> fluxsum;
            
            if (d_flux_is_face)
            {
                boost::shared_ptr<pdat::FaceDataFactory<double> > fdf(
                    BOOST_CAST<pdat::FaceDataFactory<double>,
                    hier::PatchDataFactory>(var->getPatchDataFactory()));
                TBOX_ASSERT(fdf);
                fluxsum.reset(new pdat::OuterfaceVariable<double>(
                    dim,
                    fsum_name,
                    fdf->getDepth()));
                d_flux_face_registered = true;
            }
            else
            {
                boost::shared_ptr<pdat::SideDataFactory<double> > sdf(
                    BOOST_CAST<pdat::SideDataFactory<double>,
                    hier::PatchDataFactory>(var->getPatchDataFactory()));
                TBOX_ASSERT(sdf);
                fluxsum.reset(new pdat::OutersideVariable<double>(
                    dim,
                    fsum_name,
                    sdf->getDepth()));
                d_flux_side_registered = true;
            }
            
            d_fluxsum_variables.push_back(fluxsum);
            
            int fs_id = variable_db->registerVariableAndContext(
                fluxsum,
                d_scratch,
                zero_ghosts);
            
            d_fluxsum_data.setFlag(fs_id);
            
            boost::shared_ptr<hier::CoarsenOperator> coarsen_op(
                transfer_geom->lookupCoarsenOperator(fluxsum, coarsen_name));
            
            d_coarsen_fluxsum->registerCoarsen(scr_id, fs_id, coarsen_op);
            
            break;
        }
        case SOURCE:
        {
            d_source_variables.push_back(var);
            
            int scr_id = variable_db->registerVariableAndContext(
                var,
                d_scratch,
                ghosts);
            
            d_source_var_data.setFlag(scr_id);
            
            std::vector<int> intermediate_id(d_number_steps);
            for (int sn = 0; sn < d_number_steps; sn++)
            {
                intermediate_id[sn] = variable_db->registerVariableAndContext(
                    var,
                    d_intermediate[sn],
                    ghosts_intermediate);
            }
            
            for (int sn = 0; sn < d_number_steps; sn++)
            {
                d_intermediate_source_var_data[sn].setFlag(intermediate_id[sn]);
            }
            
            break;
        }
        case TEMPORARY:
        {
            int scr_id = variable_db->registerVariableAndContext(
                var,
                d_scratch,
                ghosts);
            
            d_temp_var_scratch_data.setFlag(scr_id);
            
            // /*
            //  * Set boundary fill schedules for scratch data when data statistics are outputted.
            //  */
            
            // boost::shared_ptr<hier::RefineOperator> refine_op; // NO_REFINE
            // if (ghosts > hier::IntVector::getZero(ghosts.getDim()))
            // {
            //     d_fill_statistics->registerRefine(
            //         scr_id, scr_id, scr_id, refine_op);
            // }
            
            break;
        }
        case STATISTICS:
        {
            int scr_id = variable_db->registerVariableAndContext(
                var,
                d_scratch,
                ghosts);
            
            int stats_tmp_id = variable_db->registerVariableAndContext(
                var,
                d_stats_tmp,
                ghosts);
            
            d_stats_var_scratch_data.setFlag(scr_id);
            
            /*
             * Set boundary fill schedules for scratch data when data statistics are outputted.
             */
            
            boost::shared_ptr<hier::RefineOperator> refine_op(
                transfer_geom->lookupRefineOperator(var, refine_name));
            
            d_fill_statistics->registerRefine(
                scr_id, scr_id, stats_tmp_id, refine_op);
            
            break;
        }
        default:
        {
            TBOX_ERROR(d_object_name
                       << ":  "
                       << "unknown RK_VAR_TYPE = "
                       << RK_v_type
                       << std::endl);
        }
    }
}


/*
 **************************************************************************************************
 *
 * Process FLUX, FLUX INTEGRAL, SOURCE data before integration on the level.
 *
 * We allocate FLUX and SOURCE storage if appropriate.
 *
 * If the advance is not temporary, we also zero out the FLUX INTEGRALS on the first step of any
 * level finer than level zero.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::preprocessFluxAndSourceData(
    const boost::shared_ptr<hier::PatchLevel>& level,
    const double cur_time,
    const double new_time,
    const bool regrid_advance,
    const bool first_step,
    const bool last_step)
{
    NULL_USE(last_step);
    NULL_USE(cur_time);
    
    TBOX_ASSERT(level);
    TBOX_ASSERT(cur_time <= new_time);
    
    hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
    
    const int level_number = level->getLevelNumber();
    
    if (!regrid_advance)
    {
        if (((level_number > 0) && first_step) ||
            ((level_number == 0) && !d_have_flux_on_level_zero))
        {
            level->allocatePatchData(d_flux_var_data, new_time);
            if (level_number == 0)
            {
                d_have_flux_on_level_zero = true;
            }
            level->allocatePatchData(d_source_var_data, new_time);
        }
    }
    else
    {
        if (first_step)
        {
            level->allocatePatchData(d_flux_var_data, new_time);
            level->allocatePatchData(d_source_var_data, new_time);
        }
    }
    
    if (!regrid_advance && (level_number > 0))
    {
        if (first_step)
        {
            level->allocatePatchData(d_fluxsum_data, new_time);
            
            for (hier::PatchLevel::iterator p(level->begin());
                 p != level->end();
                 p++)
            {
                const boost::shared_ptr<hier::Patch>& patch = *p;
                
                std::list<boost::shared_ptr<hier::Variable> >::iterator fs_var =
                    d_fluxsum_variables.begin();
                
                while (fs_var != d_fluxsum_variables.end())
                {
                    int fsum_id =
                        variable_db->mapVariableAndContextToIndex(
                            *fs_var,
                            d_scratch);
                    
                    if (d_flux_is_face)
                    {
                        boost::shared_ptr<pdat::OuterfaceData<double> > fsum_data(
                            BOOST_CAST<pdat::OuterfaceData<double>, hier::PatchData>(
                                patch->getPatchData(fsum_id)));
                        
                        TBOX_ASSERT(fsum_data);
                        fsum_data->fillAll(0.0);
                    }
                    else
                    {
                        boost::shared_ptr<pdat::OutersideData<double> > fsum_data(
                            BOOST_CAST<pdat::OutersideData<double>, hier::PatchData>(
                                patch->getPatchData(fsum_id)));
                        
                        TBOX_ASSERT(fsum_data);
                        fsum_data->fillAll(0.0);
                    }
                    
                    fs_var++;
                }
            }
        }
        else
        {
            level->setTime(new_time, d_fluxsum_data);
        }
    
    } // if ( !regrid_advance && (level_number > 0) )
}


/*
 **************************************************************************************************
 *
 * Process FLUX, FLUX INTEGRAL and SOURCE data after advancing the solution on the level. During
 * normal integration steps, the flux integrals are updated for subsequent synchronization by adding
 * FLUX values to flux integrals.
 *
 * If the advance is not temporary (regular integration step):
 * 1) If the level is the finest in the hierarchy, FLUX and SOURCE data is deallocated. It is not
 *    used during synchronization, and is only maintained if needed for the advance.
 *
 * 2) If the level is not the coarsest in the hierarchy, update the  flux integrals for later
 *    synchronization by adding FLUX values to flux integrals.
 *
 * If the advance is temporary, deallocate the flux and sourc data if first step.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::postprocessFluxAndSourceData(
    const boost::shared_ptr<hier::PatchLevel>& level,
    const bool regrid_advance,
    const bool first_step,
    const bool last_step)
{
    NULL_USE(last_step);
    
    TBOX_ASSERT(level);
    
    if (level->getDim() > tbox::Dimension(3))
    {
        TBOX_ERROR("RungeKuttaLevelIntegrator::postprocessFluxAndSourceData : DIM > 3 not implemented"
            << std::endl);
    }
    
    if (regrid_advance && first_step)
    {
        level->deallocatePatchData(d_flux_var_data);
        level->deallocatePatchData(d_source_var_data);
    }
    
    if (!regrid_advance && (level->getLevelNumber() > 0))
    {
        for (hier::PatchLevel::iterator p(level->begin());
             p != level->end();
             p++)
        {
            const boost::shared_ptr<hier::Patch>& patch = *p;
            
            std::list<boost::shared_ptr<hier::Variable> >::iterator flux_var =
                d_flux_variables.begin();
            std::list<boost::shared_ptr<hier::Variable> >::iterator fluxsum_var =
                d_fluxsum_variables.begin();
            
            const hier::Index& ilo = patch->getBox().lower();
            const hier::Index& ihi = patch->getBox().upper();
            
            while (flux_var != d_flux_variables.end())
            {
                boost::shared_ptr<hier::PatchData> flux_data(
                    patch->getPatchData(*flux_var, d_scratch));
                boost::shared_ptr<hier::PatchData> fsum_data(
                   patch->getPatchData(*fluxsum_var, d_scratch));
                
                boost::shared_ptr<pdat::FaceData<double> > fflux_data;
                boost::shared_ptr<pdat::OuterfaceData<double> > ffsum_data;
                
                boost::shared_ptr<pdat::SideData<double> > sflux_data;
                boost::shared_ptr<pdat::OutersideData<double> > sfsum_data;
                
                int ddepth;
                hier::IntVector flux_num_ghosts(level->getDim());
                
                if (d_flux_is_face)
                {
                    fflux_data =
                        BOOST_CAST<pdat::FaceData<double>, hier::PatchData>(
                            flux_data);
                    ffsum_data =
                        BOOST_CAST<pdat::OuterfaceData<double>, hier::PatchData>(
                            fsum_data);
                    
                    TBOX_ASSERT(fflux_data && ffsum_data);
                    TBOX_ASSERT(fflux_data->getDepth() == ffsum_data->getDepth());
                    ddepth = fflux_data->getDepth();
                    flux_num_ghosts = fflux_data->getGhostCellWidth();
                }
                else
                {
                    sflux_data =
                        BOOST_CAST<pdat::SideData<double>, hier::PatchData>(
                            flux_data);
                    sfsum_data =
                        BOOST_CAST<pdat::OutersideData<double>, hier::PatchData>(
                            fsum_data);
                    
                    TBOX_ASSERT(sflux_data && sfsum_data);
                    TBOX_ASSERT(sflux_data->getDepth() == sfsum_data->getDepth());
                    ddepth = sflux_data->getDepth();
                    flux_num_ghosts = sflux_data->getGhostCellWidth();
                }
                
                for (int d = 0; d < ddepth; d++)
                {
                    // Loop over lower and upper parts of outer face/side arrays.
                    for (int ifs = 0; ifs < 2; ifs++)
                    {
                        if (level->getDim() == tbox::Dimension(1))
                        {
                            SAMRAI_F77_FUNC(upfluxsum1d, UPFLUXSUM1D) (
                                ilo(0),
                                ihi(0),
                                flux_num_ghosts(0),
                                ifs,
                                fflux_data->getPointer(0, d),
                                ffsum_data->getPointer(0, ifs, d));
                        }
                        else
                        {
                            if (d_flux_is_face)
                            {
                                if (level->getDim() == tbox::Dimension(2))
                                {
                                    SAMRAI_F77_FUNC(upfluxsumface2d0, UPFLUXSUMFACE2D0) (
                                        ilo(0),
                                        ilo(1),
                                        ihi(0),
                                        ihi(1),
                                        flux_num_ghosts(0),
                                        flux_num_ghosts(1),
                                        ifs,
                                        fflux_data->getPointer(0, d),
                                        ffsum_data->getPointer(0, ifs, d));
                                    SAMRAI_F77_FUNC(upfluxsumface2d1, UPFLUXSUMFACE2D1) (
                                        ilo(0),
                                        ilo(1),
                                        ihi(0),
                                        ihi(1),
                                        flux_num_ghosts(0),
                                        flux_num_ghosts(1),
                                        ifs,
                                        fflux_data->getPointer(1, d),
                                        ffsum_data->getPointer(1, ifs, d));
                                }
                                if (level->getDim() == tbox::Dimension(3))
                                {
                                    SAMRAI_F77_FUNC(upfluxsumface3d0, UPFLUXSUMFACE3D0) (
                                        ilo(0),
                                        ilo(1),
                                        ilo(2),
                                        ihi(0),
                                        ihi(1),
                                        ihi(2),
                                        flux_num_ghosts(0),
                                        flux_num_ghosts(1),
                                        flux_num_ghosts(2),
                                        ifs,
                                        fflux_data->getPointer(0, d),
                                        ffsum_data->getPointer(0, ifs, d));
                                    SAMRAI_F77_FUNC(upfluxsumface3d1, UPFLUXSUMFACE3D1) (
                                        ilo(0),
                                        ilo(1),
                                        ilo(2),
                                        ihi(0),
                                        ihi(1),
                                        ihi(2),
                                        flux_num_ghosts(0),
                                        flux_num_ghosts(1),
                                        flux_num_ghosts(2),
                                        ifs,
                                        fflux_data->getPointer(1, d),
                                        ffsum_data->getPointer(1, ifs, d));
                                    SAMRAI_F77_FUNC(upfluxsumface3d2, UPFLUXSUMFACE3D2) (
                                        ilo(0),
                                        ilo(1),
                                        ilo(2),
                                        ihi(0),
                                        ihi(1),
                                        ihi(2),
                                        flux_num_ghosts(0),
                                        flux_num_ghosts(1),
                                        flux_num_ghosts(2),
                                        ifs,
                                        fflux_data->getPointer(2, d),
                                        ffsum_data->getPointer(2, ifs, d));
                                }
                            }
                            else
                            {
                                if (level->getDim() == tbox::Dimension(2))
                                {
                                    SAMRAI_F77_FUNC(upfluxsumside2d0, UPFLUXSUMSIDE2D0) (
                                        ilo(0),
                                        ilo(1), ihi(0), ihi(1),
                                        flux_num_ghosts(0),
                                        flux_num_ghosts(1),
                                        ifs,
                                        sflux_data->getPointer(0, d),
                                        sfsum_data->getPointer(0, ifs, d));
                                    SAMRAI_F77_FUNC(upfluxsumside2d1, UPFLUXSUMSIDE2D1) (
                                        ilo(0),
                                        ilo(1), ihi(0), ihi(1),
                                        flux_num_ghosts(0),
                                        flux_num_ghosts(1),
                                        ifs,
                                        sflux_data->getPointer(1, d),
                                        sfsum_data->getPointer(1, ifs, d));
                                }
                                if (level->getDim() == tbox::Dimension(3))
                                {
                                    SAMRAI_F77_FUNC(upfluxsumside3d0, UPFLUXSUMSIDE3D0) (
                                        ilo(0),
                                        ilo(1), ilo(2),
                                        ihi(0), ihi(1), ihi(2),
                                        flux_num_ghosts(0),
                                        flux_num_ghosts(1),
                                        flux_num_ghosts(2),
                                        ifs,
                                        sflux_data->getPointer(0, d),
                                        sfsum_data->getPointer(0, ifs, d));
                                    SAMRAI_F77_FUNC(upfluxsumside3d1, UPFLUXSUMSIDE3D1) (
                                        ilo(0),
                                        ilo(1), ilo(2),
                                        ihi(0), ihi(1), ihi(2),
                                        flux_num_ghosts(0),
                                        flux_num_ghosts(1),
                                        flux_num_ghosts(2),
                                        ifs,
                                        sflux_data->getPointer(1, d),
                                        sfsum_data->getPointer(1, ifs, d));
                                    SAMRAI_F77_FUNC(upfluxsumside3d2, UPFLUXSUMSIDE3D2) (
                                        ilo(0),
                                        ilo(1), ilo(2),
                                        ihi(0), ihi(1), ihi(2),
                                        flux_num_ghosts(0),
                                        flux_num_ghosts(1),
                                        flux_num_ghosts(2),
                                        ifs,
                                        sflux_data->getPointer(2, d),
                                        sfsum_data->getPointer(2, ifs, d));
                                }
                            }  // if face operations vs. side operations
                        }
                    }  // loop over lower and upper sides/faces
                }  // loop over depth
                
                flux_var++;
                fluxsum_var++;
            }  // loop over flux variables
        }  // loop over patches
    }  // if !regrid_advance and level number > 0 ....
}


/*
 **************************************************************************************************
 *
 * Copy time-dependent data from source to destination on level.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::copyTimeDependentData(
    const boost::shared_ptr<hier::PatchLevel>& level,
    const boost::shared_ptr<hier::VariableContext>& src_context,
    const boost::shared_ptr<hier::VariableContext>& dst_context)
{
    TBOX_ASSERT(level);
    TBOX_ASSERT(src_context);
    TBOX_ASSERT(dst_context);
    
    for (hier::PatchLevel::iterator ip(level->begin());
         ip != level->end();
         ip++)
    {
        const boost::shared_ptr<hier::Patch>& patch = *ip;
        
        std::list<boost::shared_ptr<hier::Variable> >::iterator time_dep_var =
            d_time_dep_variables.begin();
            
        while (time_dep_var != d_time_dep_variables.end())
        {
            boost::shared_ptr<hier::PatchData> src_data(
                patch->getPatchData(*time_dep_var, src_context));
            
            boost::shared_ptr<hier::PatchData> dst_data(
                patch->getPatchData(*time_dep_var, dst_context));
            
            dst_data->copy(*src_data);
            
            time_dep_var++;
        }
    }
}


/*
 **************************************************************************************************
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::recordStatistics(
   const hier::PatchLevel& patch_level,
   double current_time)
{
    const int ln = patch_level.getLevelNumber();
    
    if (ln >= static_cast<int>(s_boxes_stat.size()))
    {
        s_boxes_stat.resize(ln + 1);
        s_cells_stat.resize(ln + 1);
        s_timestamp_stat.resize(ln + 1);
    }
    
    if (ln >= 0 /* Don't record work on non-hierarchy levels */)
    {
        if (!s_boxes_stat[ln])
        {
            std::string lnstr = tbox::Utilities::intToString(ln, 1);
            s_boxes_stat[ln] = tbox::Statistician::getStatistician()->
                getStatistic(std::string("HLI_BoxesL") + lnstr, "PROC_STAT");
            s_cells_stat[ln] = tbox::Statistician::getStatistician()->
                getStatistic(std::string("HLI_CellsL") + lnstr, "PROC_STAT");
            s_timestamp_stat[ln] = tbox::Statistician::getStatistician()->
                getStatistic(std::string("HLI_TimeL") + lnstr, "PROC_STAT");
        }
        
        double level_local_boxes =
            static_cast<double>(patch_level.getBoxLevel()->getLocalNumberOfBoxes());
        double level_local_cells =
            static_cast<double>(patch_level.getBoxLevel()->getLocalNumberOfCells());
        s_boxes_stat[ln]->recordProcStat(level_local_boxes);
        s_cells_stat[ln]->recordProcStat(level_local_cells);
        s_timestamp_stat[ln]->recordProcStat(current_time);
    }
}

/*
 **************************************************************************************************
 * Write out gridding statistics collected by advanceLevel.
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::printStatistics(
   std::ostream& s) const
{
    /*
     * Output statistics.
     */
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
    // Collect statistic on mesh size.
    tbox::Statistician* statn = tbox::Statistician::getStatistician();
    
    statn->finalize(false);
    // statn->printLocalStatData(s);
    if (mpi.getRank() == 0)
    {
        // statn->printAllGlobalStatData(s);
        double n_cell_updates = 0; // Number of cell updates.
        double n_patch_updates = 0; // Number of patch updates.
        for (int ln = 0; ln < static_cast<int>(s_cells_stat.size()); ln++)
        {
            tbox::Statistic& cstat = *s_cells_stat[ln];
            tbox::Statistic& bstat = *s_boxes_stat[ln];
            tbox::Statistic& tstat = *s_timestamp_stat[ln];
            s << "statistic " << cstat.getName() << ":" << std::endl;
            if (0)
            {
                s << "Global: \n";
                statn->printGlobalProcStatDataFormatted(cstat.getInstanceId(), s);
            }
            s << "Seq#   SimTime           C-Sum   C-Avg   C-Min ->      C-Max  C-Max/Avg     B-Sum    "
              << "B-Avg B-Min -> B-Max B-Max/Avg  C/B-Avg\n";
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
            for (int sn = 0; sn < cstat.getStatSequenceLength(); sn++)
            {
                const double csum = statn->getGlobalProcStatSum(
                    cstat.getInstanceId(),
                    sn);
                const double cmax = statn->getGlobalProcStatMax(
                    cstat.getInstanceId(),
                    sn);
                const double cmin = statn->getGlobalProcStatMin(
                    cstat.getInstanceId(),
                    sn);
                const double cavg = csum / mpi.getSize();
                const double cmaxnorm = cavg != 0 ? cmax / cavg : 1;
                const double bsum = statn->getGlobalProcStatSum(
                    bstat.getInstanceId(),
                    sn);
                const double bmax = statn->getGlobalProcStatMax(
                    bstat.getInstanceId(),
                    sn);
                const double bmin = statn->getGlobalProcStatMin(
                    bstat.getInstanceId(),
                    sn);
                const double bavg = bsum / mpi.getSize();
                const double bmaxnorm = bavg != 0 ? bmax / bavg : 1;
                const double stime = statn->getGlobalProcStatMin(
                    tstat.getInstanceId(),
                    sn);
                s << std::setw(4) << sn
                  << " " << std::scientific << std::setprecision(6) << std::setw(12) << stime
                  << " " << std::fixed << std::setprecision(0) << std::setw(12) << csum
                  << " " << std::setw(7) << cavg
                  << " " << std::setw(7) << cmin
                  << " -> " << std::setw(10) << cmax
                  << "  " << std::setw(9) << std::setprecision(2) << cmaxnorm
                  << " " << std::fixed << std::setprecision(0)
                  << std::setw(9) << bsum
                  << "  " << std::fixed << std::setprecision(2)
                  << std::setw(7) << bavg
                  << "   " << std::fixed << std::setprecision(0)
                  << std::setw(3) << bmin
                  << " -> " << std::setw(5) << bmax
                  << "  " << std::setw(8) << std::setprecision(2) << bmaxnorm
                  << "   " << std::setw(6) << std::setprecision(0) << (bsum != 0 ? csum / bsum : 0)
                  << std::endl;
                n_cell_updates += csum;
                n_patch_updates += bsum;
            }
        }
        s << "Total number of cell updates: "
          << n_cell_updates
          << " (" << double(n_cell_updates) / mpi.getSize()
          << "/proc)"
          << std::endl;
        s << "Total number of boxe updates: "
          << n_patch_updates
          << " (" << double(n_patch_updates) / mpi.getSize()
          << "/proc)"
          << std::endl;
    }
}


/*
 **************************************************************************************************
 *
 * Print all class data for RungeKuttaLevelIntegrator object.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::printClassData(
   std::ostream& os) const
{
    os << "\nPrint RungeKuttaLevelIntegrator object..." << std::endl;
    os << std::endl;
    
    os << "RungeKuttaLevelIntegrator: this = "
       << (RungeKuttaLevelIntegrator *)this << std::endl;
    os << "d_object_name = " << d_object_name << std::endl;
    
    os << "d_use_cfl = " << d_use_cfl << std::endl;
    
    if (d_use_cfl)
    {
        os << "d_cfl = " << d_cfl << "\n"
           << "d_cfl_init = " << d_cfl_init << std::endl;
        os << "d_lag_dt_computation = " << d_lag_dt_computation << "\n"
           << "d_use_ghosts_for_dt = "
           << d_use_ghosts_for_dt << std::endl;
    }
    else
    {
        os << "d_dt = " << d_dt << std::endl;
    }
    
    os << "d_patch_strategy = "
       << (RungeKuttaPatchStrategy *)d_patch_strategy << std::endl;
    
    os << std::endl;
    
    os << "Runge-Kutta weights:" << std::endl;
    for (int i = 0; i < static_cast<int>(d_alpha.size()); i++)
    {
        if (i == 0)
        {
            os << "   d_alpha = ";
        }
        else
        {
            os << "             ";
        }
        
        os << d_alpha[i][0];
        for (int j = 1; j < static_cast<int>(d_alpha[i].size()); j++)
        {
            os << " , " << d_alpha[i][j];
        }
        os << std::endl;
    }
    for (int i = 0; i < static_cast<int>(d_beta.size()); i++)
    {
        if (i == 0)
        {
            os << "   d_beta  = ";
        }
        else
        {
            os << "             ";
        }
        
        os << d_beta[i][0];
        for (int j = 1; j < static_cast<int>(d_beta[i].size()); j++)
        {
            os << " , " << d_beta[i][j];
        }
        os << std::endl;
    }
    for (int i = 0; i < static_cast<int>(d_gamma.size()); i++)
    {
        if (i == 0)
        {
            os << "   d_gamma = ";
        }
        else
        {
            os << "             ";
        }
        
        os << d_gamma[i][0];
        for (int j = 1; j < static_cast<int>(d_gamma[i].size()); j++)
        {
            os << " , " << d_gamma[i][j];
        }
        os << std::endl;
    }
    os << std::endl;
    os << "NOTE: Not printing variable arrays, ComponentSelectors, communication schedules, etc."
       << std::endl;
}


/*
 **************************************************************************************************
 *
 * Writes out the class version number, d_cfl, d_cfl_init, d_lag_dt_computation, and d_use_ghosts_for_dt
 * to the restart database.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::putToRestart(
    const boost::shared_ptr<tbox::Database>& restart_db) const
{
    TBOX_ASSERT(restart_db);
    
    restart_db->putInteger("RUNGE_KUTTA_LEVEL_INTEGRATOR_VERSION",
        RUNGE_KUTTA_LEVEL_INTEGRATOR_VERSION);
    
    restart_db->putBool("use_cfl", d_use_cfl);
    restart_db->putDouble("cfl", d_cfl);
    restart_db->putDouble("cfl_init", d_cfl_init);
    restart_db->putBool("lag_dt_computation", d_lag_dt_computation);
    restart_db->putBool("use_ghosts_to_compute_dt", d_use_ghosts_for_dt);
    restart_db->putDouble("dt", d_dt);
    restart_db->putBool("DEV_distinguish_mpi_reduction_costs",
        d_distinguish_mpi_reduction_costs);
    
    restart_db->putDatabase("RungeKuttaWeights");
    
    boost::shared_ptr<tbox::Database> RK_db(restart_db->getDatabase("RungeKuttaWeights"));
    
    RK_db->putInteger("number_steps", d_number_steps);
    
    for (int i = 0; i < d_number_steps; i++)
    {
        std::string alpha_array_name = "alpha_" + tbox::Utilities::intToString(i);
        
        std::vector<double> alpha_array(static_cast<int>(d_alpha[i].size()));
        
        for (int j = 0; j < static_cast<int>(d_alpha[i].size()); j++)
        {
            alpha_array[j] = d_alpha[i][j];
        }
        
        RK_db->putDoubleVector(alpha_array_name, alpha_array);
    }
    
    for (int i = 0; i < d_number_steps; i++)
    {
        std::string beta_array_name = "beta_" + tbox::Utilities::intToString(i);
        
        std::vector<double> beta_array(static_cast<int>(d_beta[i].size()));
        
        for (int j = 0; j < static_cast<int>(d_beta[i].size()); j++)
        {
            beta_array[j] = d_beta[i][j];
        }
        
        RK_db->putDoubleVector(beta_array_name, beta_array);
    }
    
    for (int i = 0; i < d_number_steps; i++)
    {
        std::string gamma_array_name = "gamma_" + tbox::Utilities::intToString(i);
        
        std::vector<double> gamma_array(static_cast<int>(d_gamma[i].size()));
        
        for (int j = 0; j < static_cast<int>(d_gamma[i].size()); j++)
        {
            gamma_array[j] = d_gamma[i][j];
        }
        
        RK_db->putDoubleVector(gamma_array_name, gamma_array);
    }
}


/*
 **************************************************************************************************
 *
 * Output statistics of data.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::outputDataStatistics(
    const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const double statistics_data_time)
{
    t_output_data_statistics->start();
    
    int num_levels = hierarchy->getNumberOfLevels();
    
    for (int li = 0; li < num_levels; li++)
    {
        boost::shared_ptr<hier::PatchLevel> patch_level(
            hierarchy->getPatchLevel(li));
        
        patch_level->allocatePatchData(d_saved_var_scratch_data, statistics_data_time);
        // patch_level->allocatePatchData(d_temp_var_scratch_data, statistics_data_time);
        patch_level->allocatePatchData(d_stats_var_scratch_data, statistics_data_time);
    }
    
    d_patch_strategy->setDataContext(d_scratch);
    
    for (int li = 0; li < num_levels; li++)
    {
        d_bdry_sched_advance[li]->fillData(statistics_data_time);
    }
    
    // Compute variables if necessary.
    
    d_patch_strategy->computeStatisticsVariables(hierarchy);
    
    // Exchange halo values of variables if necessary.
    
    for (int li = 0; li < num_levels; li++)
    {
        boost::shared_ptr<hier::PatchLevel> patch_level(
            hierarchy->getPatchLevel(li));
        
        boost::shared_ptr<xfer::RefineSchedule> fill_schedule;
        
        if (li > 0)
        {
            fill_schedule = d_fill_statistics->createSchedule(
                patch_level,
                li - 1,
                hierarchy);
        }
        else
        {
            fill_schedule = d_fill_statistics->createSchedule(
                patch_level,
                d_patch_strategy);
        }
        
        fill_schedule->fillData(statistics_data_time);
    }
    
    // Filter variables if necessary.
    
    const int num_filtering = 4;
    
    for (int count = 0; count < num_filtering; count++)
    {
        d_patch_strategy->filterStatisticsVariables(hierarchy);
        
        // Exchange halo values of variables if necessary.
        
        for (int li = 0; li < num_levels; li++)
        {
            boost::shared_ptr<hier::PatchLevel> patch_level(
                hierarchy->getPatchLevel(li));
            
            boost::shared_ptr<xfer::RefineSchedule> fill_schedule;
            
            if (li > 0)
            {
                fill_schedule = d_fill_statistics->createSchedule(
                    patch_level,
                    li - 1,
                    hierarchy);
            }
            else
            {
                fill_schedule = d_fill_statistics->createSchedule(
                    patch_level,
                    d_patch_strategy);
            }
            
            fill_schedule->fillData(statistics_data_time);
        }
    }
    
    // Set conservative variables as the filtered valued if necessary.
    // Then, compute and output statistics.
    // Reset conservative variables as the pre-filtered values if necessary.
    
    d_patch_strategy->outputDataStatistics(hierarchy, statistics_data_time);
    
    for (int li = 0; li < num_levels; li++)
    {
        boost::shared_ptr<hier::PatchLevel> patch_level(
            hierarchy->getPatchLevel(li));
        
        patch_level->deallocatePatchData(d_stats_var_scratch_data);
        // patch_level->deallocatePatchData(d_temp_var_scratch_data);
        patch_level->deallocatePatchData(d_saved_var_scratch_data);
    }
    
    t_output_data_statistics->stop();
}


/*
 **************************************************************************************************
 *
 * Reads in cfl, cfl_init, lag_dt_computation, and use_ghosts_to_compute_dt from the input database.
 * Note all restart values are overriden with values from the input database.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::getFromInput(
    const boost::shared_ptr<tbox::Database>& input_db,
    bool is_from_restart)
{
    if (!is_from_restart && !input_db)
    {
        TBOX_ERROR(": RungeKuttaLevelIntegrator::getFromInput()\n"
                   << "no input database supplied" << std::endl);
    }
    
    if (!is_from_restart)
    {
        d_use_cfl = input_db->getBoolWithDefault("use_cfl", true);
        
        if (d_use_cfl)
        {
            d_cfl = input_db->getDouble("cfl");
            d_cfl_init = input_db->getDouble("cfl_init");
            d_lag_dt_computation = input_db->getBoolWithDefault("lag_dt_computation", true);
            d_use_ghosts_for_dt = input_db->getBoolWithDefault("use_ghosts_to_compute_dt", false);
        }
        else
        {
            d_dt = input_db->getDouble("dt");
        }
        
        d_distinguish_mpi_reduction_costs = input_db->getBoolWithDefault("DEV_distinguish_mpi_reduction_costs", false);
        
        if (input_db->keyExists("RungeKuttaWeights"))
        {
            boost::shared_ptr<tbox::Database> RK_db(input_db->getDatabase("RungeKuttaWeights"));
            
            d_number_steps = RK_db->getInteger("number_steps");
            
            d_alpha.resize(d_number_steps);
            d_beta.resize(d_number_steps);
            d_gamma.resize(d_number_steps);
            
            for (int sn = 0; sn < d_number_steps; sn++)
            {
                std::string alpha_array_name = "alpha_" + tbox::Utilities::intToString(sn);
                std::string beta_array_name = "beta_" + tbox::Utilities::intToString(sn);
                std::string gamma_array_name = "gamma_" + tbox::Utilities::intToString(sn);
                
                if (RK_db->keyExists(alpha_array_name))
                {
                    size_t alpha_array_size = RK_db->getArraySize(alpha_array_name);
                    
                    if (static_cast<int>(alpha_array_size) == sn + 1)
                    {
                        d_alpha[sn] = RK_db->getDoubleVector(alpha_array_name);
                    }
                    else
                    {
                        TBOX_ERROR(": RungeKuttaLevelIntegrator::getFromInput()\n"
                                   << "number of '"
                                   << alpha_array_name
                                   << "' entries must be = "
                                   << (sn + 1)
                                   << std::endl);
                    }
                }
                else
                {
                    TBOX_ERROR(": RungeKuttaLevelIntegrator::getFromInput()\n"
                               << "Key data '"
                               << alpha_array_name
                               << "' not found in input. "
                               << std::endl);
                }
                
                if (RK_db->keyExists(beta_array_name))
                {
                    size_t beta_array_size = RK_db->getArraySize(beta_array_name);
                    
                    if (static_cast<int>(beta_array_size) == sn + 1)
                    {
                        d_beta[sn] = RK_db->getDoubleVector(beta_array_name);
                    }
                    else
                    {
                        TBOX_ERROR(": RungeKuttaLevelIntegrator::getFromInput()\n"
                                   << "number of '"
                                   << beta_array_name
                                   << "' entries must be = "
                                   << (sn + 1)
                                   << std::endl);
                    }
                }
                else
                {
                    TBOX_ERROR(": RungeKuttaLevelIntegrator::getFromInput()\n"
                               << "Key data '"
                               << beta_array_name
                               << "' not found in input. "
                               << std::endl);
                }
                
                if (RK_db->keyExists(gamma_array_name))
                {
                    size_t gamma_array_size = RK_db->getArraySize(gamma_array_name);
                    
                    if (static_cast<int>(gamma_array_size) == sn + 1)
                    {
                        d_gamma[sn] = RK_db->getDoubleVector(gamma_array_name);
                    }
                    else
                    {
                        TBOX_ERROR(": RungeKuttaLevelIntegrator::getFromInput()\n"
                                   << "number of '"
                                   << gamma_array_name
                                   << "' entries must be = "
                                   << (sn + 1)
                                   << std::endl);
                    }
                }
                else
                {
                    TBOX_ERROR(": RungeKuttaLevelIntegrator::getFromInput()\n"
                               << "Key data '"
                               << gamma_array_name
                               << "' not found in input. "
                               << std::endl);
                }
            }
        }
        else    // Use SSPRK(3, 3) Runge-Kutta scheme as the default scheme.
        {
            d_number_steps = 3;
            
            d_alpha.resize(d_number_steps);
            d_beta.resize(d_number_steps);
            d_gamma.resize(d_number_steps);
            
            for (int sn = 0; sn < d_number_steps; sn++)
            {
                d_alpha[sn].resize(sn + 1);
                d_beta[sn].resize(sn + 1);
                d_gamma[sn].resize(sn + 1);
            }
            
            d_alpha[0][0] = 1.0;
            d_alpha[1][0] = 3.0/4.0;
            d_alpha[1][1] = 1.0/4.0;
            d_alpha[2][0] = 1.0/3.0;
            d_alpha[2][1] = 0.0;
            d_alpha[2][2] = 2.0/3.0;
            
            d_beta[0][0]  = 1.0;
            d_beta[1][0]  = 0.0;
            d_beta[1][1]  = 1.0/4.0;
            d_beta[2][0]  = 0.0;
            d_beta[2][1]  = 0.0;
            d_beta[2][2]  = 2.0/3.0;
            
            d_gamma[0][0] = 1.0/6.0;
            d_gamma[1][0] = 0.0;
            d_gamma[1][1] = 1.0/6.0;
            d_gamma[2][0] = 0.0;
            d_gamma[2][1] = 0.0;
            d_gamma[2][2] = 2.0/3.0;
        }
    }
    else if (input_db)
    {
        bool read_on_restart = input_db->getBoolWithDefault("read_on_restart", false);
     
        if (read_on_restart)
        {
            d_use_cfl = input_db->getBoolWithDefault("use_cfl", true);
            
            if (d_use_cfl)
            {
                d_cfl = input_db->getDoubleWithDefault("cfl", d_cfl);
                
                d_cfl_init = input_db->getDoubleWithDefault("cfl_init", d_cfl_init);
                
                d_lag_dt_computation =
                   input_db->getBoolWithDefault("lag_dt_computation",
                      d_lag_dt_computation);
                
                d_use_ghosts_for_dt =
                   input_db->getBoolWithDefault("use_ghosts_to_compute_dt",
                      d_use_ghosts_for_dt);
            }
            else
            {
                d_dt = input_db->getDouble("dt");
            }
            
            d_distinguish_mpi_reduction_costs =
                input_db->getBoolWithDefault("DEV_distinguish_mpi_reduction_costs",
                    d_distinguish_mpi_reduction_costs);
        }
    }
}


/*
 **************************************************************************************************
 *
 * First, gets the database corresponding to the object_name from the restart file. If this database
 * exists, this method checks to make sure that the version number of the class matches the version
 * number of the restart file. If they match, then d_cfl, d_cfl_init, d_lag_dt_computation, and
 * d_use_ghosts_to_compute_dt are read from restart database. Note all restart values can be overriden
 * with values from the input database.
 *
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::getFromRestart()
{
    boost::shared_ptr<tbox::Database> root_db(
        tbox::RestartManager::getManager()->getRootDatabase());

    if (!root_db->isDatabase(d_object_name))
    {
        TBOX_ERROR("Restart database corresponding to "
                   << d_object_name
                   << " not found in restart file" << std::endl);
    }
    boost::shared_ptr<tbox::Database> db(root_db->getDatabase(d_object_name));
    
    int ver = db->getInteger("RUNGE_KUTTA_LEVEL_INTEGRATOR_VERSION");
    if (ver != RUNGE_KUTTA_LEVEL_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name
                   << ":  "
                   << "Restart file version different "
                   << "than class version."
                   << std::endl);
    }
    
    d_use_cfl = db->getBool("use_cfl");
    d_cfl = db->getDouble("cfl");
    d_cfl_init = db->getDouble("cfl_init");
    d_lag_dt_computation = db->getBool("lag_dt_computation");
    d_use_ghosts_for_dt = db->getBool("use_ghosts_to_compute_dt");
    d_dt = db->getDouble("dt");
    d_distinguish_mpi_reduction_costs = db->getBool("DEV_distinguish_mpi_reduction_costs");
    
    boost::shared_ptr<tbox::Database> RK_db(db->getDatabase("RungeKuttaWeights"));
    
    d_number_steps = RK_db->getInteger("number_steps");
    
    d_alpha.resize(d_number_steps);
    d_beta.resize(d_number_steps);
    d_gamma.resize(d_number_steps);
    
    for (int sn = 0; sn < d_number_steps; sn++)
    {
        std::string alpha_array_name = "alpha_" + tbox::Utilities::intToString(sn);
        std::string beta_array_name = "beta_" + tbox::Utilities::intToString(sn);
        std::string gamma_array_name = "gamma_" + tbox::Utilities::intToString(sn);
        
        d_alpha[sn] = RK_db->getDoubleVector(alpha_array_name);
        d_beta[sn] = RK_db->getDoubleVector(beta_array_name);
        d_gamma[sn] = RK_db->getDoubleVector(gamma_array_name);
    }
}


/*
 **************************************************************************************************
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::initializeCallback()
{
    /*
     * Timers:  one for each of the communication algorithms ("create"
     * indicates schedule creation, "comm" indicates communication)
     */
    t_advance_bdry_fill_comm = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::advance_bdry_fill_comm");
    t_error_bdry_fill_create = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::error_bdry_fill_create");
    t_error_bdry_fill_comm = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::error_bdry_fill_comm");
    t_barrier_after_error_bdry_fill_comm = tbox::TimerManager::getManager()->
        getTimer( "RungeKuttaLevelIntegrator::barrier_after_error_bdry_fill_comm");
    t_advance_mpi_reductions = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::advance_mpi_reductions");
    t_get_level_dt_mpi_reductions = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::get_level_dt_mpi_reductions");
    t_initialize_level_data = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::initializeLevelData()");
    t_init_level_create_sched = tbox::TimerManager::getManager()->
        getTimer( "RungeKuttaLevelIntegrator::initializeLevelData()_createSched");
    t_init_level_fill_data = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::initializeLevelData()_fillData");
    t_init_level_fill_interior = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::initializeLevelData()_fill_interior");
    t_advance_bdry_fill_create = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::advance_bdry_fill_create");
    t_new_advance_bdry_fill_create = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::new_advance_bdry_fill_create");
    t_apply_value_detector = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::applyValueDetector()");
    t_apply_gradient_detector = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::applyGradientDetector()");
    t_apply_multiresolution_detector = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::applyMultiresolutionDetector()");
    t_apply_integral_detector = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::applyIntegralDetector()");
    t_tag_cells = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::tag_cells");
    t_coarsen_rich_extrap = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::coarsen_rich_extrap");
    t_get_level_dt = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::getLevelDt()");
    t_get_level_dt_sync = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::getLevelDt()_sync");
    t_advance_level = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::advanceLevel()");
    t_new_advance_bdry_fill_comm = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::new_advance_bdry_fill_comm");
    t_patch_num_kernel = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::patch_numerical_kernels");
    t_advance_level_sync = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::advanceLevel()_sync");
    t_std_level_sync = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::standardLevelSynchronization()");
    t_sync_new_levels = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::synchronizeNewLevels()");
    t_sync_initial_create = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::sync_initial_create");
    t_sync_initial_comm = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::sync_initial_comm");
    t_coarsen_fluxsum_create = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::coarsen_fluxsum_create");
    t_coarsen_fluxsum_comm = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::coarsen_fluxsum_comm");
    t_coarsen_sync_create = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::coarsen_sync_create");
    t_coarsen_sync_comm = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::coarsen_sync_comm");
    t_output_data_statistics = tbox::TimerManager::getManager()->
        getTimer("RungeKuttaLevelIntegrator::output_data_statistics");
}


/*
 **************************************************************************************************
 **************************************************************************************************
 */
void
RungeKuttaLevelIntegrator::finalizeCallback()
{
    t_advance_bdry_fill_comm.reset();
    t_error_bdry_fill_create.reset();
    t_error_bdry_fill_comm.reset();
    t_barrier_after_error_bdry_fill_comm.reset();
    t_advance_mpi_reductions.reset();
    t_get_level_dt_mpi_reductions.reset();
    t_initialize_level_data.reset();
    t_init_level_create_sched.reset();
    t_init_level_fill_data.reset();
    t_init_level_fill_interior.reset();
    t_advance_bdry_fill_create.reset();
    t_new_advance_bdry_fill_create.reset();
    t_apply_value_detector.reset();
    t_apply_gradient_detector.reset();
    t_apply_multiresolution_detector.reset();
    t_apply_integral_detector.reset();
    t_tag_cells.reset();
    t_coarsen_rich_extrap.reset();
    t_get_level_dt.reset();
    t_get_level_dt_sync.reset();
    t_advance_level.reset();
    t_new_advance_bdry_fill_comm.reset();
    t_patch_num_kernel.reset();
    t_advance_level_sync.reset();
    t_std_level_sync.reset();
    t_sync_new_levels.reset();
    t_sync_initial_create.reset();
    t_sync_initial_comm.reset();
    t_coarsen_fluxsum_create.reset();
    t_coarsen_fluxsum_comm.reset();
    t_coarsen_sync_create.reset();
    t_coarsen_sync_comm.reset();
    t_output_data_statistics.reset();

#ifdef HLI_RECORD_STATS
    /*
     * Statistics on number of cells and patches generated.
     */
    s_boxes_stat.clear();
    s_cells_stat.clear();
    s_timestamp_stat.clear();
#endif
}

