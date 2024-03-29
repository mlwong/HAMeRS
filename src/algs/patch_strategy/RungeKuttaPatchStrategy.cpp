/************************************************************************
 *
 * This file is modified from HyperbolicPatchStrategy.c of the SAMRAI
 * distribution. For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Interface to patch routines for hyperbolic Runge-Kutta
 *                integration scheme.
 *
 ************************************************************************/

#include "algs/patch_strategy/RungeKuttaPatchStrategy.hpp"

#include "SAMRAI/tbox/Utilities.h"

RungeKuttaPatchStrategy::RungeKuttaPatchStrategy():
    xfer::RefinePatchStrategy(),
    xfer::CoarsenPatchStrategy(),
    d_data_context()
{
}


RungeKuttaPatchStrategy::~RungeKuttaPatchStrategy()
{
}


/*
 **************************************************************************************************
 *
 * Default virtual function implementations.
 *
 **************************************************************************************************
 */

void
RungeKuttaPatchStrategy::preprocessTagCellsValueDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const double regrid_time,
    const bool initial_error,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(level_number);
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::postprocessTagCellsValueDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const double regrid_time,
    const bool initial_error,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(level_number);
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::tagCellsOnPatchValueDetector(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_index,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(patch);
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(tag_index);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    TBOX_ERROR("RungeKuttaPatchStrategy::tagCellsOnPatchValueDetector()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method."
        << std::endl);
}


void
RungeKuttaPatchStrategy::preprocessTagCellsGradientDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const double regrid_time,
    const bool initial_error,
    const bool uses_value_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(level_number);
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::postprocessTagCellsGradientDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const double regrid_time,
    const bool initial_error,
    const bool uses_value_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(level_number);
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::tagCellsOnPatchGradientDetector(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_index,
    const bool uses_value_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(patch);
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(tag_index);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    TBOX_ERROR("RungeKuttaPatchStrategy::tagCellsOnPatchGradientDetector()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method."
        << std::endl);
}


void
RungeKuttaPatchStrategy::preprocessTagCellsMultiresolutionDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const double regrid_time,
    const bool initial_error,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(level_number);
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::postprocessTagCellsMultiresolutionDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const double regrid_time,
    const bool initial_error,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(level_number);
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::tagCellsOnPatchMultiresolutionDetector(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_index,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(patch);
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(tag_index);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    TBOX_ERROR("RungeKuttaPatchStrategy::tagCellsOnPatchMultiresolutionDetector()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method."
        << std::endl);
}


void
RungeKuttaPatchStrategy::preprocessTagCellsIntegralDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const double regrid_time,
    const bool initial_error,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(level_number);
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::postprocessTagCellsIntegralDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const int level_number,
    const double regrid_time,
    const bool initial_error,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(level_number);
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::tagCellsOnPatchIntegralDetector(
    hier::Patch& patch,
    const double regrid_time,
    const bool initial_error,
    const int tag_index,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(patch);
    NULL_USE(regrid_time);
    NULL_USE(initial_error);
    NULL_USE(tag_index);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    TBOX_ERROR("RungeKuttaPatchStrategy::tagCellsOnPatchIntegralDetector()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method."
        << std::endl);
}


void
RungeKuttaPatchStrategy::tagCellsOnPatchRichardsonExtrapolation(
    hier::Patch& patch,
    const int error_level_number,
    const HAMERS_SHARED_PTR<hier::VariableContext>& coarsened_fine,
    const HAMERS_SHARED_PTR<hier::VariableContext>& advanced_coarse,
    const double regrid_time,
    const double deltat,
    const int error_coarsen_ratio,
    const bool initial_error,
    const int tag_index,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too)
{
    NULL_USE(patch);
    NULL_USE(error_level_number);
    NULL_USE(coarsened_fine);
    NULL_USE(advanced_coarse);
    NULL_USE(regrid_time);
    NULL_USE(deltat);
    NULL_USE(error_coarsen_ratio);
    NULL_USE(initial_error);
    NULL_USE(tag_index);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    TBOX_ERROR("RungeKuttaPatchStrategy::tagRichardsonExtrapolationCells()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method." << std::endl);
}


/**
 * Compute variables for computing the statistics of data.
 */
void
RungeKuttaPatchStrategy::computeStatisticsVariables(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy)
{
    NULL_USE(patch_hierarchy);
}


/**
 * Filter variables for computing the statistics of data.
 */
void
RungeKuttaPatchStrategy::filterStatisticsVariables(
    const int level,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy)
{
    NULL_USE(level);
    NULL_USE(patch_hierarchy);
}


/**
 * Output the header of statistics.
 */
void
RungeKuttaPatchStrategy::outputHeaderStatistics()
{
}


/**
 * Compute the statistics of data.
 */
void
RungeKuttaPatchStrategy::computeDataStatistics(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double statistics_data_time)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(statistics_data_time);
}


/**
 * Output the statistics of data.
 */
void
RungeKuttaPatchStrategy::outputDataStatistics(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
    const double output_time)
{
    NULL_USE(patch_hierarchy);
    NULL_USE(output_time);
}


void
RungeKuttaPatchStrategy::setupLoadBalancer(
    RungeKuttaLevelIntegrator* integrator,
    mesh::GriddingAlgorithm* gridding_algorithm)
{
    NULL_USE(integrator);
    NULL_USE(gridding_algorithm);
}


void
RungeKuttaPatchStrategy::preprocessAdvanceLevelState(
    const HAMERS_SHARED_PTR<hier::PatchLevel>& level,
    double current_time,
    double dt,
    bool first_step,
    bool last_step,
    bool regrid_advance)
{
    NULL_USE(level);
    NULL_USE(current_time);
    NULL_USE(dt);
    NULL_USE(first_step);
    NULL_USE(last_step);
    NULL_USE(regrid_advance);
}


void
RungeKuttaPatchStrategy::postprocessAdvanceLevelState(
    const HAMERS_SHARED_PTR<hier::PatchLevel>& level,
    double current_time,
    double dt,
    bool first_step,
    bool last_step,
    bool regrid_advance)
{
    NULL_USE(level);
    NULL_USE(current_time);
    NULL_USE(dt);
    NULL_USE(first_step);
    NULL_USE(last_step);
    NULL_USE(regrid_advance);
}
