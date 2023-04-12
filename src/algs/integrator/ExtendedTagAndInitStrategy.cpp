/*************************************************************************
 *
 * This file is modified from StandardTagAndInitStrategy.c of the SAMRAI
 * distribution. For full copyright information, see COPYRIGHT and
 * COPYING.LESSER of SAMRAI distribution.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for error detection.
 *
 ************************************************************************/

#include "algs/integrator/ExtendedTagAndInitStrategy.hpp"

#include "SAMRAI/tbox/Utilities.h"

ExtendedTagAndInitStrategy::ExtendedTagAndInitStrategy()
{
}


ExtendedTagAndInitStrategy::~ExtendedTagAndInitStrategy()
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
ExtendedTagAndInitStrategy::applyRefineRegions(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(hierarchy);
    NULL_USE(level_number);
    NULL_USE(error_data_time);
    NULL_USE(tag_index);
    NULL_USE(initial_time);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    TBOX_ERROR("ExtendedTagAndInitStrategy::applyRefineRegions()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method." << std::endl);
}


/*
 *************************************************************************
 *
 * Default virtual function implementations.
 *
 *************************************************************************
 */

void
ExtendedTagAndInitStrategy::applyValueDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_refine_regions_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(hierarchy);
    NULL_USE(level_number);
    NULL_USE(error_data_time);
    NULL_USE(tag_index);
    NULL_USE(initial_time);
    NULL_USE(uses_refine_regions_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    TBOX_ERROR("ExtendedTagAndInitStrategy::applyValueDetector()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method." << std::endl);
}


void
ExtendedTagAndInitStrategy::applyGradientDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_refine_regions_too,
    const bool uses_value_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(hierarchy);
    NULL_USE(level_number);
    NULL_USE(error_data_time);
    NULL_USE(tag_index);
    NULL_USE(initial_time);
    NULL_USE(uses_refine_regions_too);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    TBOX_ERROR("ExtendedTagAndInitStrategy::applyGradientDetector()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method." << std::endl);
}


void
ExtendedTagAndInitStrategy::applyMultiresolutionDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_refine_regions_too,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_integral_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(hierarchy);
    NULL_USE(level_number);
    NULL_USE(error_data_time);
    NULL_USE(tag_index);
    NULL_USE(initial_time);
    NULL_USE(uses_refine_regions_too);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_integral_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    TBOX_ERROR("ExtendedTagAndInitStrategy::applyMultiresolutionDetector()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method." << std::endl);
}


void
ExtendedTagAndInitStrategy::applyIntegralDetector(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_refine_regions_too,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_richardson_extrapolation_too)
{
    NULL_USE(hierarchy);
    NULL_USE(level_number);
    NULL_USE(error_data_time);
    NULL_USE(tag_index);
    NULL_USE(initial_time);
    NULL_USE(uses_refine_regions_too);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_richardson_extrapolation_too);
    TBOX_ERROR("ExtendedTagAndInitStrategy::applyIntegralDetector()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method." << std::endl);
}


void
ExtendedTagAndInitStrategy::applyRichardsonExtrapolation(
    const HAMERS_SHARED_PTR<hier::PatchLevel>& level,
    const double error_data_time,
    const int tag_index,
    const double deltat,
    const int error_coarsen_ratio,
    const bool initial_time,
    const bool uses_refine_regions_too,
    const bool uses_value_detector_too,
    const bool uses_gradient_detector_too,
    const bool uses_multiresolution_detector_too,
    const bool uses_integral_detector_too)
{
    NULL_USE(level);
    NULL_USE(error_data_time);
    NULL_USE(tag_index);
    NULL_USE(deltat);
    NULL_USE(error_coarsen_ratio);
    NULL_USE(initial_time);
    NULL_USE(uses_refine_regions_too);
    NULL_USE(uses_value_detector_too);
    NULL_USE(uses_gradient_detector_too);
    NULL_USE(uses_multiresolution_detector_too);
    NULL_USE(uses_integral_detector_too);
    TBOX_ERROR("ExtendedTagAndInitStrategy::applyRichardsonExtrapolation()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method." << std::endl);
}


void
ExtendedTagAndInitStrategy::coarsenDataForRichardsonExtrapolation(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const HAMERS_SHARED_PTR<hier::PatchLevel>& coarser_level,
    const double coarsen_data_time,
    const bool before_advance)
{
    NULL_USE(hierarchy);
    NULL_USE(level_number);
    NULL_USE(coarser_level);
    NULL_USE(coarsen_data_time);
    NULL_USE(before_advance);
    TBOX_ERROR("ExtendedTagAndInitStrategy::coarsenDataForRichardsonExtrapolation()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method." << std::endl);
}


double
ExtendedTagAndInitStrategy::getLevelDt(
    const HAMERS_SHARED_PTR<hier::PatchLevel>& level,
    const double dt_time,
    const bool initial_time)
{
    NULL_USE(level);
    NULL_USE(dt_time);
    NULL_USE(initial_time);
    TBOX_ERROR("ExtendedTagAndInitStrategy::getLevelDt()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method." << std::endl);
    
    return 0.0;
}


void
ExtendedTagAndInitStrategy::resetTimeDependentData(
   const HAMERS_SHARED_PTR<hier::PatchLevel>& level,
   const double new_time,
   const bool can_be_refined)
{
   NULL_USE(level);
   NULL_USE(new_time);
   NULL_USE(can_be_refined);
   TBOX_ERROR("ExtendedTagAndInitStrategy::resetTimeDependentData()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method." << std::endl);
}


double
ExtendedTagAndInitStrategy::advanceLevel(
    const HAMERS_SHARED_PTR<hier::PatchLevel>& level,
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
    const double current_time,
    const double new_time,
    const bool first_step,
    const bool last_step,
    const bool regrid_advance)
{
    NULL_USE(level);
    NULL_USE(hierarchy);
    NULL_USE(current_time);
    NULL_USE(new_time);
    NULL_USE(first_step);
    NULL_USE(last_step);
    NULL_USE(regrid_advance);
    TBOX_ERROR("ExtendedTagAndInitStrategy::advanceLevel()"
        << "\nNo derived class supplies a concrete implementation for "
        << "\nthis method." << std::endl);
    
    return 0.0;
}


void
ExtendedTagAndInitStrategy::resetDataToPreadvanceState(
    const HAMERS_SHARED_PTR<hier::PatchLevel>& level)
{
    NULL_USE(level);
    TBOX_ERROR("ExtendedTagAndInitStrategy::resetDataToPreadvanceState()"
        << "\nNo class derived supplies a concrete implementation for "
        << "\nthis method." << std::endl);
}


void
ExtendedTagAndInitStrategy::processHierarchyBeforeAddingNewLevel(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const HAMERS_SHARED_PTR<hier::BoxLevel>& new_box_level)
{
    NULL_USE(hierarchy);
    NULL_USE(level_number);
    NULL_USE(new_box_level);
}


void
ExtendedTagAndInitStrategy::processLevelBeforeRemoval(
    const HAMERS_SHARED_PTR<hier::PatchHierarchy>& hierarchy,
    const int level_number,
    const HAMERS_SHARED_PTR<hier::PatchLevel>& old_level)
{
    NULL_USE(hierarchy);
    NULL_USE(level_number);
    NULL_USE(old_level);
}

