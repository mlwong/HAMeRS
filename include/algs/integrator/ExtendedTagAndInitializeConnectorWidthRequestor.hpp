/*************************************************************************
 *
 * This file is modified from StandardTagAndInitializeConnectorWidthRequestor.h
 * of the SAMRAI version 3.9.1 distribution.  For full copyright information,
 * see COPYRIGHT and COPYING.LESSER of the SAMRAI distribution.
 *
 ************************************************************************/

#ifndef EXTENDED_TAG_AND_INITIALIZE_CONNECTOR_WIDTH_REQUESTOR_HPP
#define EXTENDED_TAG_AND_INITIALIZE_CONNECTOR_WIDTH_REQUESTOR_HPP

#include "HAMeRS_config.hpp"

#include "SAMRAI/hier/PatchHierarchy.h"

using namespace SAMRAI;

/*!
 * @brief Implementation of the strategy class
 * hier::PatchHierarchy::ConnectorWidthRequestorStrategy to tell the
 * hier::PatchHierarchy how wide ExtendedTagAndInitialize needs
 * Connectors between hierarchy levels to be.
 *
 * To do Richardson extrapolation, ExtendedTagAndInitialize will
 * coarsen a level and populate it with data.  A coarsened level has a
 * bigger ghost region because the coarse cells are bigger.  This
 * class is for telling the PatchHierarchy that
 * ExtendedTagAndInitialize will request Connectors based on the width
 * of the coarsened level.
 */

class ExtendedTagAndInitializeConnectorWidthRequestor:
    public hier::PatchHierarchy::ConnectorWidthRequestorStrategy
{
    public:
    /*!
     * @brief Constructor.
     */
    ExtendedTagAndInitializeConnectorWidthRequestor();
    
    /*!
     * @brief Compute Connector widths that this class requires in
     * order to work properly on a given hierarchy.
     *
     * Implements the pure virtual method
     * hier::PatchHierarchy::ConnectorWidthRequestorStrategy::computeRequiredConnectorWidths()
     *
     * @param[out] self_connector_widths Array of widths for Connectors
     * from a level to itself.
     *
     * @param[out] fine_connector_widths Array of widths for Connectors
     * from a level to the next finer level.
     *
     * @param[in]  patch_hierarchy
     */
    void
    computeRequiredConnectorWidths(
        std::vector<hier::IntVector>& self_connector_widths,
        std::vector<hier::IntVector>& fine_connector_widths,
        const hier::PatchHierarchy& patch_hierarchy) const;
    
    private:
    /*!
     * @brief Return the coarsen ratio to be used with Richardson
     * extrapolation on the PatchHierarchy.
     *
     * @param[in] ratios_to_coarser Refinement ratios in a hierarchy.
     * @c ratios_to_coarser[ln] is the ratio between level ln and level
     * ln-1.
     *
     * @pre (ratios_to_coarser[1](0) % 2 == 0) ||
     *      (ratios_to_coarser[1](0) % 3 == 0)
     */
    int
    computeCoarsenRatio(
        const std::vector<hier::IntVector>& ratios_to_coarser) const;
    
};

#endif /* EXTENDED_TAG_AND_INITIALIZE_CONNECTOR_WIDTH_REQUESTOR_HPP */
