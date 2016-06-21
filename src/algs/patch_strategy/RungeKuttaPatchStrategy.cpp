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
 *************************************************************************
 *
 * Default virtual function implementations.
 *
 *************************************************************************
 */

void
RungeKuttaPatchStrategy::tagGradientDetectorCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_index,
   const bool uses_multiresolution_detector_too,
   const bool uses_integral_detector_too,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(patch);
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(tag_index);
   NULL_USE(uses_multiresolution_detector_too);
   NULL_USE(uses_integral_detector_too);
   NULL_USE(uses_richardson_extrapolation_too);
   TBOX_ERROR("RungeKuttaPatchStrategy::tagGradientDetectorCells()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method."
      << std::endl);
}


void
RungeKuttaPatchStrategy::preprocessTagGradientDetectorCells(
   const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
   const int level_number,
   const double regrid_time,
   const bool initial_error,
   const bool uses_multiresolution_detector_too,
   const bool uses_integral_detector_too,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(patch_hierarchy);
   NULL_USE(level_number);
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(uses_multiresolution_detector_too);
   NULL_USE(uses_integral_detector_too);
   NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::postprocessTagGradientDetectorCells(
   const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
   const int level_number,
   const double regrid_time,
   const bool initial_error,
   const bool uses_multiresolution_detector_too,
   const bool uses_integral_detector_too,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(patch_hierarchy);
   NULL_USE(level_number);
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(uses_multiresolution_detector_too);
   NULL_USE(uses_integral_detector_too);
   NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::tagMultiresolutionDetectorCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_index,
   const bool uses_gradient_detector_too,
   const bool uses_integral_detector_too,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(patch);
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(tag_index);
   NULL_USE(uses_gradient_detector_too);
   NULL_USE(uses_integral_detector_too);
   NULL_USE(uses_richardson_extrapolation_too);
   TBOX_ERROR("RungeKuttaPatchStrategy::tagMultiresolutionDetectorCells()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method."
      << std::endl);
}


void
RungeKuttaPatchStrategy::preprocessTagMultiresolutionDetectorCells(
   const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
   const int level_number,
   const double regrid_time,
   const bool initial_error,
   const bool uses_gradient_detector_too,
   const bool uses_integral_detector_too,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(patch_hierarchy);
   NULL_USE(level_number);
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(uses_gradient_detector_too);
   NULL_USE(uses_integral_detector_too);
   NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::postprocessTagMultiresolutionDetectorCells(
   const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
   const int level_number,
   const double regrid_time,
   const bool initial_error,
   const bool uses_gradient_detector_too,
   const bool uses_integral_detector_too,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(patch_hierarchy);
   NULL_USE(level_number);
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(uses_gradient_detector_too);
   NULL_USE(uses_integral_detector_too);
   NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::tagIntegralDetectorCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_index,
   const bool uses_gradient_detector_too,
   const bool uses_multiresolution_detector_too,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(patch);
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(tag_index);
   NULL_USE(uses_gradient_detector_too);
   NULL_USE(uses_multiresolution_detector_too);
   NULL_USE(uses_richardson_extrapolation_too);
   TBOX_ERROR("RungeKuttaPatchStrategy::tagIntegralDetectorCells()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method."
      << std::endl);
}


void
RungeKuttaPatchStrategy::preprocessTagIntegralDetectorCells(
   const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
   const int level_number,
   const double regrid_time,
   const bool initial_error,
   const bool uses_gradient_detector_too,
   const bool uses_multiresolution_detector_too,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(patch_hierarchy);
   NULL_USE(level_number);
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(uses_gradient_detector_too);
   NULL_USE(uses_multiresolution_detector_too);
   NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::postprocessTagIntegralDetectorCells(
   const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
   const int level_number,
   const double regrid_time,
   const bool initial_error,
   const bool uses_gradient_detector_too,
   const bool uses_multiresolution_detector_too,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(patch_hierarchy);
   NULL_USE(level_number);
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(uses_gradient_detector_too);
   NULL_USE(uses_multiresolution_detector_too);
   NULL_USE(uses_richardson_extrapolation_too);
}


void
RungeKuttaPatchStrategy::tagRichardsonExtrapolationCells(
   hier::Patch& patch,
   const int error_level_number,
   const boost::shared_ptr<hier::VariableContext>& coarsened_fine,
   const boost::shared_ptr<hier::VariableContext>& advanced_coarse,
   const double regrid_time,
   const double deltat,
   const int error_coarsen_ratio,
   const bool initial_error,
   const int tag_index,
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
   NULL_USE(uses_gradient_detector_too);
   NULL_USE(uses_multiresolution_detector_too);
   NULL_USE(uses_integral_detector_too);
   TBOX_ERROR("RungeKuttaPatchStrategy::tagRichardsonExtrapolationCells()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method." << std::endl);
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
   const boost::shared_ptr<hier::PatchLevel>& level,
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
   const boost::shared_ptr<hier::PatchLevel>& level,
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
