/*************************************************************************
 *
 * This file is modified from StandardTagAndInitializeConnectorWidthRequestor.c
 * of the SAMRAI distribution.  For full copyright information, see COPYRIGHT
 * and COPYING.LESSER of SAMRAI distribution.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   ExtendedTagAndInitialize's implementation of PatchHierarchy
 *
 ************************************************************************/

#include "algs/integrator/ExtendedTagAndInitializeConnectorWidthRequestor.hpp"

#include "SAMRAI/xfer/RefineScheduleConnectorWidthRequestor.h"
#include "SAMRAI/tbox/MathUtilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

/*
 *************************************************************************
 * Static function that computes greatest common divisor.
 *************************************************************************
 */
static int GCD(
    const int a,
    const int b)
{
    int at = tbox::MathUtilities<int>::Min(a, b);
    int bt = tbox::MathUtilities<int>::Max(a, b);
    
    if (at == 0 || bt == 0) return bt;
    
    at = (at > 0 ? at : -at);
    bt = (bt > 0 ? bt : -bt);
    
    int r0 = at;
    int r1 = at;
    int r2 = bt;
    
    while (!(r2 == 0))
    {
        r0 = r1;
        r1 = r2;
        int q = r0 / r1;
        r2 = r0 - r1 * q;
    }
    
    return r1;
}


/*
 **************************************************************************
 **************************************************************************
 */
ExtendedTagAndInitializeConnectorWidthRequestor::ExtendedTagAndInitializeConnectorWidthRequestor()
{
}


/*
 **************************************************************************
 * Compute Connector widths that this class requires in order to work
 * properly on a given hierarchy.
 *
 * ExtendedTagAndInitialize generates coarsened versions of levels on
 * the hierarchy and uses RefineSchedule to fill its ghost data.
 * It is the RefineSchedule that drives the need for Connector width
 * when using ExtendedTagAndInitialize.
 *
 * Coarsening in Richardson extrapolation has the effect of making the
 * ghost cells look bigger when RefineSchedule computes data
 * dependency.  The required Connector widths are equivalent to those
 * that RefineSchedule generates for the additional ghost widths.
 * Therefore, this class requires the same Connector widths that
 * RefineSchedule requires, except for bigger ghost cells.
 *
 **************************************************************************
 */
void
ExtendedTagAndInitializeConnectorWidthRequestor::computeRequiredConnectorWidths(
    std::vector<hier::IntVector>& self_connector_widths,
    std::vector<hier::IntVector>& fine_connector_widths,
    const hier::PatchHierarchy& patch_hierarchy) const
{
    const tbox::Dimension& dim(patch_hierarchy.getDim());
    
    /*
     * Get the refinement ratios on the hierarchy.  This is the ratio
     * by which ExtendedTagAndInitialize may coarsen a level in the
     * hierarchy.  It is the growth factor for ghost data.
     */
    std::vector<hier::IntVector> ratios_to_coarser(
        patch_hierarchy.getMaxNumberOfLevels(),
        hier::IntVector(dim));
    
    for (int ln = 0; ln < patch_hierarchy.getMaxNumberOfLevels(); ln++)
    {
        ratios_to_coarser[ln] = patch_hierarchy.getRatioToCoarserLevel(ln);
    }
    
    int error_coarsen_ratio = computeCoarsenRatio(ratios_to_coarser);
    
    /*
     * Use the error coarsen ratio to compute the Connector widths
     * associated with the coarsened ghost regions.
     */
    xfer::RefineScheduleConnectorWidthRequestor rscwri;
    rscwri.setGhostCellWidthFactor(error_coarsen_ratio);
    rscwri.computeRequiredConnectorWidths(
        self_connector_widths,
        fine_connector_widths,
        patch_hierarchy);
}


/*
 *************************************************************************
 * Compute Error coarsen ratio for Richardson extrapolation. For a given
 * level, the error coarsen ratio should be the greatest common divisor
 * (GCD) of the refinement ratio applied to the level.  This value
 * should generally be 2 or 3 (e.g. refinement ratio=2 gives GCD=2;
 * rr=3 gives GCD=3; rr=4 gives GCD=2; etc.).
 *
 * Note that this algorithm was lifted from ExtendedTagAndInitialize.
 *************************************************************************
 */
int
ExtendedTagAndInitializeConnectorWidthRequestor::computeCoarsenRatio(
    const std::vector<hier::IntVector>& ratios_to_coarser) const
{
    const tbox::Dimension& dim(ratios_to_coarser[0].getDim());
    /*
     * Compute GCD on first coordinate direction of level 1
     */
    int error_coarsen_ratio = 0;
    int gcd_level1 = ratios_to_coarser[1](0);
    if ((gcd_level1 % 2) == 0)
    {
        error_coarsen_ratio = 2;
    }
    else if ((gcd_level1 % 3) == 0)
    {
        error_coarsen_ratio = 3;
    }
    else
    {
        TBOX_ERROR("Unable to perform Richardson extrapolation algorithm "
            << "with ratios_to_coarser[1](0) = " << gcd_level1 << "\n"
            << "Did you intend to use Richardson extrapolation?\n"
            << "If no, you don't have to be here.  If yes, fix\n"
            << "the refinement ratios in your hierarchy.");
    }
    
    /*
     * Iterate through levels and check the coarsen ratios to make sure the
     * error coarsen ratios computed in every coordinate direction on every
     * level are between the supported 2 or 3, and that the error coarsen
     * ratios are constant over the hierarchy.
     */
    for (int ln = 1; ln < static_cast<int>(ratios_to_coarser.size()); ln++)
    {
        for (int d = 0; d < dim.getValue(); d++)
        {
            int gcd = GCD(error_coarsen_ratio, ratios_to_coarser[ln](d));
            if ((gcd % error_coarsen_ratio) != 0)
            {
                gcd = ratios_to_coarser[ln](d);
                TBOX_ERROR("ExtendedTagAndInitializeConnectorWidthRequestor::computeCoarsenRatio:\n"
                    << "Unable to perform Richardson extrapolation because\n"
                    << "the error coarsen ratio computed from the\n"
                    << "ratios_to_coarser entries is not constant across all\n"
                    << "levels, in all coordinate directions, of the hierarchy. In\n"
                    << "order to use Richardson extrapolation, the minimum\n"
                    << "divisor (> 1) of all the ratios_to_coarser entries must\n"
                    << "be 2 -or- 3:\n"
                    << "   level 1(0): minimum divisor: "
                    << error_coarsen_ratio
                    << "\n   level " << ln << "(" << d
                    << "):"
                    << ": ratios_to_coarser = " << gcd << "\n"
                    << "Did you intend to use Richardson extrapolation?\n"
                    << "If no, you don't have to be here.  If yes, fix\n"
                    << "the refinement ratios in your hierarchy.");
            }
        }
    }
    
    return error_coarsen_ratio;
}
