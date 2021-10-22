#include "util/basic_boundary_conditions/BasicCartesianBoundaryUtilities1.hpp"

#include "util/basic_boundary_conditions/BasicBoundaryConditions.hpp"
#include "util/basic_boundary_conditions/CartesianBoundaryDefines.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <algorithm>

// Integer constant for debugging improperly set boundary data.
#define BOGUS_BDRY_LOC (-9999)

/*
 * This function reads 1D boundary data from given input database.
 * The integer boundary condition types are placed in the integer
 * arrays supplied by the caller (typically, the concrete BoundaryStrategy
 * provided).  When DIRICHLET or NEUMANN conditions are specified, control
 * is passed to the BoundaryStrategy to read the boundary state data
 * specific to the problem.
 *
 * Errors will be reported and the program will abort whenever necessary
 * boundary condition information is missing in the input database, or
 * when the data read in is either unknown or inconsistent.  The periodic
 * domain information is used to determine which boundary node entries are
 * not required from input.
 *
 * Arguments are:
 *    bdry_strategy .... object that reads DIRICHLET or NEUMANN data
 *    input_db ......... input database containing all boundary data
 *    node_locs ........ array of locations of nodes for applying
 *                       boundary conditions.
 *    node_conds ....... array into which integer boundary conditions
 *                       for nodes are read
 *    periodic ......... integer vector specifying which coordinate
 *                       directions are periodic (value returned from
 *                       GridGeometry2::getPeriodicShift())
 */
void
BasicCartesianBoundaryUtilities1::getFromInput(
    BoundaryUtilityStrategy* bdry_strategy,
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    std::vector<int>& node_locs,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(1));
    TBOX_ASSERT(bdry_strategy != 0);
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_1D_NODES);
    if (static_cast<int>(node_locs.size()) > 0)
    {
        TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
        TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_1D_NODES);
    }
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_1D_NODES);
    
    if (!input_db)
    {
        TBOX_ERROR("BasicCartesianBoundaryUtilities1::getFromInput()\n"
            << "No input database supplied."
            << std::endl);
    }
    
    if (static_cast<int>(node_locs.size()) > 0)
    {
        read1dBdryNodes(
            bdry_strategy,
            input_db,
            node_locs,
            node_conds,
            periodic);
    }
}


/*
 * Function to remove 1d boundary nodes with boundary conditions filled by this class for a patch.
 *
 * Arguments are:
 *    bdry_node_locs ....... array of locations of nodes for applying
 *                           boundary conditions.
 *    patch ................ patch on which data object lives
 *    bdry_node_conds ...... array of boundary conditions for patch nodes
 */
void
BasicCartesianBoundaryUtilities1::removeBoundaryNodeLocations(
    std::vector<int>& bdry_node_locs,
    const hier::Patch& patch,
    const std::vector<int>& bdry_node_conds)
{
    TBOX_ASSERT(static_cast<int>(bdry_node_locs.size()) <= NUM_1D_NODES);
    TBOX_ASSERT(*min_element(bdry_node_locs.begin(), bdry_node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_node_locs.begin(), bdry_node_locs.end()) < NUM_1D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_1D_NODES);
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    const std::vector<hier::BoundaryBox>& node_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::NODE1D);
    
    for (int ni = 0; ni < static_cast<int>(node_bdry.size()); ni++)
    {
        TBOX_ASSERT(node_bdry[ni].getBoundaryType() == BDRY::NODE1D);
        
        int node_loc = node_bdry[ni].getLocationIndex();
        
        if (std::find(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc) !=
            bdry_node_locs.end())
        {
            if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::DIRICHLET)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::NEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::FLOW)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::REFLECT)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::SYMMETRY)
            {
                // Remove node locations that have boundary conditions identified.
                bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                    bdry_node_locs.end());
            }
        }
    }
}


/*
 * Function to fill node boundary values.
 *
 * Arguments are:
 *    var_name ............. name of variable (for error reporting)
 *    var_data ............. cell-centered patch data object to check
 *    patch ................ patch on which data object lives
 *    bdry_node_locs ....... array of locations of nodes for applying
 *                           boundary conditions.
 *    bdry_node_conds ...... array of boundary conditions for patch nodes
 *    bdry_node_values ..... array of boundary values for nodes
 *                           (this must be consistent with boundary
 *                           condition types)
 *    ghost_width_to_fill .. width of ghost region to fill
 */
void
BasicCartesianBoundaryUtilities1::fillNodeBoundaryData(
    const std::string& var_name,
    const HAMERS_SHARED_PTR<pdat::CellData<double> >& var_data,
    const hier::Patch& patch,
    const std::vector<int>& bdry_node_locs,
    const std::vector<int>& bdry_node_conds,
    const std::vector<double>& bdry_node_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(!var_name.empty());
    TBOX_ASSERT(var_data);
    TBOX_ASSERT(static_cast<int>(bdry_node_locs.size()) <= NUM_1D_NODES);
    TBOX_ASSERT(*min_element(bdry_node_locs.begin(), bdry_node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_node_locs.begin(), bdry_node_locs.end()) < NUM_1D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_1D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_node_values.size()) == NUM_1D_NODES*(var_data->getDepth()));
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(1));
    TBOX_ASSERT_OBJDIM_EQUALITY3(*var_data, patch, ghost_width_to_fill);
    
    NULL_USE(var_name);
    
    const HAMERS_SHARED_PTR<geom::CartesianPatchGeometry> patch_geom(
        HAMERS_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    const hier::IntVector& num_ghosts(var_data->getGhostCellWidth());
    
    /*
     * Determine the ghost cell width to fill.
     */
    hier::IntVector gcw_to_fill(tbox::Dimension(1));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(1)))
    {
        gcw_to_fill = var_data->getGhostCellWidth();
    }
    else
    {
        gcw_to_fill = hier::IntVector::min(
            num_ghosts,
            ghost_width_to_fill);
    }
    
    // Get the dimensions of box that covers the interior of patch.
    const hier::Box& interior_box(patch.getBox());
    hier::Index interior_box_lo_idx(interior_box.lower());
    hier::Index interior_box_hi_idx(interior_box.upper());
    
    /*
     * Offset the indices.
     */
    interior_box_lo_idx = interior_box_lo_idx - interior_box.lower();
    interior_box_hi_idx = interior_box_hi_idx - interior_box.lower();
    
    // Get the dimensions of box that covers interior of patch plus
    // ghost cells.
    const hier::IntVector ghostcell_dims = var_data->getGhostBox().numberCells();
    
    const std::vector<hier::BoundaryBox>& node_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::NODE1D);
    
    const int var_depth = var_data->getDepth();
    
    for (int ni = 0; ni < static_cast<int>(node_bdry.size()); ni++)
    {
        TBOX_ASSERT(node_bdry[ni].getBoundaryType() == BDRY::NODE1D);
        
        int node_loc = node_bdry[ni].getLocationIndex();
        
        if (std::find(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc) !=
            bdry_node_locs.end())
        {
            hier::Box fill_box(patch_geom->getBoundaryFillBox(
                node_bdry[ni],
                interior_box,
                gcw_to_fill));
            
            hier::Index fill_box_lo_idx(fill_box.lower());
            hier::Index fill_box_hi_idx(fill_box.upper());
            
            /*
             * Offset the indices.
             */
            fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
            fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
            
            if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::DIRICHLET)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    const int idx_cell = i + num_ghosts[0];
                    
                    for (int di = 0; di < var_depth; di++)
                    {
                        var_data->getPointer(di)[idx_cell] = bdry_node_values[node_loc*var_depth + di];
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::NEUMANN)
            {
                // NOT YET IMPLEMENTED
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::FLOW)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    const int idx_cell = i + num_ghosts[0];
                    
                    int idx_cell_pivot = idx_cell;
                    
                    if (node_loc == BDRY_LOC::XLO)
                    {
                        idx_cell_pivot = interior_box_lo_idx[0] + num_ghosts[0];
                    }
                    else if (node_loc == BDRY_LOC::XHI)
                    {
                        idx_cell_pivot = interior_box_hi_idx[0] + num_ghosts[0];
                    }
                    
                    for (int di = 0; di < var_depth; di++)
                    {
                        var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                    }
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::REFLECT)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    const int idx_cell = i + num_ghosts[0];
                    
                    int idx_cell_pivot = idx_cell;
                    
                    if (node_loc == BDRY_LOC::XLO)
                    {
                        idx_cell_pivot = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0];
                    }
                    else if (node_loc == BDRY_LOC::XHI)
                    {
                        idx_cell_pivot = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0];
                    }
                    
                    for (int di = 0; di < var_depth; di++)
                    {
                        var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                    }
                    
                    var_data->getPointer(0)[idx_cell] = -var_data->getPointer(0)[idx_cell_pivot];
                }
            }
            else if (bdry_node_conds[node_loc] == BDRY_COND::BASIC::SYMMETRY)
            {
                for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                {
                    const int idx_cell = i + num_ghosts[0];
                    
                    int idx_cell_pivot = idx_cell;
                    
                    if (node_loc == BDRY_LOC::XLO)
                    {
                        idx_cell_pivot = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) + num_ghosts[0];
                    }
                    else if (node_loc == BDRY_LOC::XHI)
                    {
                        idx_cell_pivot = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) + num_ghosts[0];
                    }
                    
                    for (int di = 0; di < var_depth; di++)
                    {
                        var_data->getPointer(di)[idx_cell] = var_data->getPointer(di)[idx_cell_pivot];
                    }
                }
            }
        }
    }
}


/*
 * Private function to read 1D node boundary data from input database.
 */
void
BasicCartesianBoundaryUtilities1::read1dBdryNodes(
    BoundaryUtilityStrategy* bdry_strategy,
    const HAMERS_SHARED_PTR<tbox::Database>& input_db,
    std::vector<int>& node_locs,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(1));
    
    TBOX_ASSERT(bdry_strategy != 0);
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_1D_NODES);
    TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_1D_NODES);
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_1D_NODES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 1; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 1)
    {
        // node boundary input required
        for (int ni = 0; ni < static_cast<int>(node_locs.size()); ni++)
        {
            int s = node_locs[ni];
            
            std::string bdry_loc_str;
            switch (s)
            {
                case BDRY_LOC::XLO:
                {
                    bdry_loc_str = "boundary_node_xlo";
                    break;
                }
                case BDRY_LOC::XHI:
                {
                    bdry_loc_str = "boundary_node_xhi";
                    break;
                }
                default: NULL_STATEMENT;
            }
            
            bool need_data_read = true;
            if (num_per_dirs > 0)
            {
                if (periodic(0) && (s == BDRY_LOC::XLO || s == BDRY_LOC::XHI))
                {
                    need_data_read = false;
                }
            }
            
            if (need_data_read)
            {
                HAMERS_SHARED_PTR<tbox::Database> bdry_loc_db(
                    input_db->getDatabase(bdry_loc_str));
                std::string bdry_cond_str =
                    bdry_loc_db->getString("boundary_condition");
                if (bdry_cond_str == "FLOW")
                {
                    node_conds[s] = BDRY_COND::BASIC::FLOW;
                    
                    node_locs[ni] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "REFLECT")
                {
                    node_conds[s] = BDRY_COND::BASIC::REFLECT;
                    
                    node_locs[ni] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "SYMMETRY")
                {
                    node_conds[s] = BDRY_COND::BASIC::SYMMETRY;
                    
                    node_locs[ni] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "DIRICHLET")
                {
                    node_conds[s] = BDRY_COND::BASIC::DIRICHLET;
                    bdry_strategy->readDirichletBoundaryDataEntry(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    node_locs[ni] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "NEUMANN")
                {
                    node_conds[s] = BDRY_COND::BASIC::NEUMANN;
                    bdry_strategy->readNeumannBoundaryDataEntry(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    node_locs[ni] = BOGUS_BDRY_LOC;
                }
            } // if (need_data_read)
       } // for (int ni = 0 ...
    } // if (num_per_dirs < 1)
    
    // Remove node locations that have boundary conditions identified.
    node_locs.erase(std::remove(node_locs.begin(), node_locs.end(), BOGUS_BDRY_LOC), node_locs.end());
}
