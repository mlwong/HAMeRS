#include "flow/flow_models/single-species/FlowModelBoundaryUtilitiesSingleSpecies.hpp"

#include "util/basic_boundary_conditions/CartesianBoundaryDefines.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

// Integer constant for debugging improperly set boundary data.
#define BOGUS_BDRY_LOC (-9999)

FlowModelBoundaryUtilitiesSingleSpecies::FlowModelBoundaryUtilitiesSingleSpecies(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const int& num_species,
    const int& num_eqn,
    const boost::shared_ptr<EquationOfStateMixingRules>& equation_of_state_mixing_rules):
        FlowModelBoundaryUtilities(
            object_name,
            dim,
            num_species,
            num_eqn,
            equation_of_state_mixing_rules)
{
    std::vector<double*> thermo_properties_ptr;
    
    const int num_thermo_properties = d_equation_of_state_mixing_rules->
        getNumberOfSpeciesThermodynamicProperties();
    
    thermo_properties_ptr.reserve(num_thermo_properties);
    d_thermo_properties.resize(num_thermo_properties);
    
    for (int ti = 0; ti < num_thermo_properties; ti++)
    {
        thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
    }
    
    d_equation_of_state_mixing_rules->getSpeciesThermodynamicProperties(
        thermo_properties_ptr,
        0);
    
    /*
     * Set default sizes for containers of boundary conditions.
     */
    
    if (d_dim == tbox::Dimension(1))
    {
        d_bdry_node_adiabatic_no_slip_vel.resize(NUM_1D_NODES);
        
        d_bdry_node_isothermal_no_slip_T.resize(NUM_1D_NODES);
        d_bdry_node_isothermal_no_slip_vel.resize(NUM_1D_NODES);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        d_bdry_edge_adiabatic_no_slip_vel.resize(NUM_2D_EDGES*2);
        
        d_bdry_edge_isothermal_no_slip_T.resize(NUM_2D_EDGES);
        d_bdry_edge_isothermal_no_slip_vel.resize(NUM_2D_EDGES*2);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        d_bdry_face_adiabatic_no_slip_vel.resize(NUM_3D_FACES*3);
        
        d_bdry_face_isothermal_no_slip_T.resize(NUM_3D_FACES);
        d_bdry_face_isothermal_no_slip_vel.resize(NUM_3D_FACES*3);
    }
}


/*
 * Function to read 1d boundary data from input database.
 * Node locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::getFromInput1d(
    const boost::shared_ptr<tbox::Database>& input_db,
    std::vector<int>& node_locs,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(1));
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_1D_NODES);
    TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_1D_NODES);
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_1D_NODES);
    
    if (!input_db)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::getFromInput1d()\n"
            << "no input database supplied"
            << std::endl);
    }
    
    read1dBdryNodes(
        input_db,
        node_locs,
        node_conds,
        periodic);
}


/*
 * Function to read 2d boundary data from input database.
 * Node and edge locations that have boundary conditions identified are removed from the
 * containers.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::getFromInput2d(
    const boost::shared_ptr<tbox::Database>& input_db,
    std::vector<int>& edge_locs,
    std::vector<int>& node_locs,
    std::vector<int>& edge_conds,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));
    TBOX_ASSERT(static_cast<int>(edge_locs.size()) <= NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_2D_NODES);
    TBOX_ASSERT(*min_element(edge_locs.begin(), edge_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(edge_locs.begin(), edge_locs.end()) < NUM_2D_EDGES);
    TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_2D_NODES);
    TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_2D_NODES);
    
    if (!input_db)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::getFromInput2d()\n"
            << "no input database supplied"
            << std::endl);
    }
    
    read2dBdryEdges(
        input_db,
        edge_locs,
        edge_conds,
        periodic);
    
    read2dBdryNodes(
        input_db,
        node_locs,
        edge_conds,
        node_conds,
        periodic);
}


/*
 * Function to read 3d boundary data from input database.
 * Node, edge and face locations that have boundary conditions identified are removed from
 * the containers.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::getFromInput3d(
    const boost::shared_ptr<tbox::Database>& input_db,
    std::vector<int>& face_locs,
    std::vector<int>& edge_locs,
    std::vector<int>& node_locs,
    std::vector<int>& face_conds,
    std::vector<int>& edge_conds,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(3));
    TBOX_ASSERT(static_cast<int>(face_locs.size()) <= NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(edge_locs.size()) <= NUM_3D_EDGES);
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_3D_NODES);
    TBOX_ASSERT(*min_element(face_locs.begin(), face_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(face_locs.begin(), face_locs.end()) < NUM_3D_FACES);
    TBOX_ASSERT(*min_element(edge_locs.begin(), edge_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(edge_locs.begin(), edge_locs.end()) < NUM_3D_EDGES);
    TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_3D_NODES);
    TBOX_ASSERT(static_cast<int>(face_conds.size()) == NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_3D_EDGES);
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_3D_NODES);
    
    if (!input_db)
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::getFromInput3d()\n"
            << "no input database supplied"
            << std::endl);
    }
    
    read3dBdryFaces(
        input_db,
        face_locs,
        face_conds,
        periodic);
    
    read3dBdryEdges(
        input_db,
        edge_locs,
        face_conds,
        edge_conds,
        periodic);
    
    read3dBdryNodes(
        input_db,
        node_locs,
        face_conds,
        node_conds,
        periodic);
}


/*
 * Function to fill 1d node boundary values for a patch.
 * Node locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::fill1dNodeBoundaryData(
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
    const hier::Patch& patch,
    std::vector<int>& bdry_node_locs,
    const std::vector<int>& bdry_node_conds,
    const std::vector<std::vector<double> >& bdry_node_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(static_cast<int>(conservative_var_data.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT(conservative_var_data[vi]);
    }
    TBOX_ASSERT(static_cast<int>(bdry_node_locs.size()) <= NUM_1D_NODES);
    TBOX_ASSERT(*min_element(bdry_node_locs.begin(), bdry_node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_node_locs.begin(), bdry_node_locs.end()) < NUM_1D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_1D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_node_values.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(bdry_node_values.size()); vi++)
    {
        TBOX_ASSERT(static_cast<int>(bdry_node_values[vi].size()) ==
                    NUM_1D_NODES*(conservative_var_data[vi]->getDepth()));
    }
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(1));
    
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY3(*conservative_var_data[vi], patch, ghost_width_to_fill);
    }
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    hier::IntVector num_ghosts(conservative_var_data[0]->getGhostCellWidth());
    for (int vi = 1; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        num_ghosts = hier::IntVector::min(
            num_ghosts,
            conservative_var_data[vi]->getGhostCellWidth());
    }
    
    /*
     * Determine the ghost cell width to fill.
     */
    
    hier::IntVector gcw_to_fill(tbox::Dimension(1));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(1)))
    {
        gcw_to_fill = num_ghosts;
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
    
    const std::vector<hier::BoundaryBox>& node_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::NODE1D);
    
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
            
            if ((bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP))
            {
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<double*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const double*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        const int idx_cell_rho = i + num_subghosts_conservative_var[0][0];
                        const int idx_cell_mom = i + num_subghosts_conservative_var[1][0];
                        const int idx_cell_E = i + num_subghosts_conservative_var[2][0];
                        
                        int idx_cell_pivot_rho = idx_cell_rho;
                        int idx_cell_pivot_mom = idx_cell_mom;
                        int idx_cell_pivot_E = idx_cell_E;
                        
                        if (node_loc == BDRY_LOC::XLO)
                        {
                            idx_cell_pivot_rho = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[0][0];
                            
                            idx_cell_pivot_mom = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[1][0];
                            
                            idx_cell_pivot_E = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[2][0];
                            
                        }
                        else if (node_loc == BDRY_LOC::XHI)
                        {
                            idx_cell_pivot_rho = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[0][0];
                            
                            idx_cell_pivot_mom = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[1][0];
                            
                            idx_cell_pivot_E = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[2][0];
                        }
                        
                        /*
                         * Set the values for density and momentum.
                         */
                        
                        Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                        Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                            2.0*Q[0][idx_cell_pivot_rho]*d_bdry_node_adiabatic_no_slip_vel[node_loc];
                        
                        /*
                         * Set the values for total internal energy.
                         */
                        
                        double epsilon_pivot = (Q[2][idx_cell_pivot_E] -
                            0.5*Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom]/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                        
                        double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &Q[0][idx_cell_pivot_rho],
                                &epsilon_pivot,
                                thermo_properties_ptr);
                        
                        double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getTemperature(
                                &Q[0][idx_cell_pivot_rho],
                                &p_pivot,
                                thermo_properties_ptr);
                        
                        double T = T_pivot;
                        
                        double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getInternalEnergyFromTemperature(
                                &Q[0][idx_cell_rho],
                                &T,
                                thermo_properties_ptr);
                        
                        double E = Q[0][idx_cell_rho]*epsilon +
                            0.5*Q[1][idx_cell_mom]*Q[1][idx_cell_mom]/Q[0][idx_cell_rho];
                        
                        Q[2][idx_cell_E] = E;
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        const int idx_cell_rho = i + num_subghosts_conservative_var[0][0];
                        const int idx_cell_mom = i + num_subghosts_conservative_var[1][0];
                        const int idx_cell_E   = i + num_subghosts_conservative_var[2][0];
                        
                        int idx_cell_pivot_rho = idx_cell_rho;
                        int idx_cell_pivot_mom = idx_cell_mom;
                        int idx_cell_pivot_E = idx_cell_E;
                        
                        if (node_loc == BDRY_LOC::XLO)
                        {
                            idx_cell_pivot_rho = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[0][0];
                            
                            idx_cell_pivot_mom = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[1][0];
                            
                            idx_cell_pivot_E = interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                num_subghosts_conservative_var[2][0];
                        }
                        else if (node_loc == BDRY_LOC::XHI)
                        {
                            idx_cell_pivot_rho = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[0][0];
                            
                            idx_cell_pivot_mom = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[1][0];
                            
                            idx_cell_pivot_E = interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                num_subghosts_conservative_var[2][0];
                        }
                        
                        /*
                         * Set the values for density, momentum and total internal energy.
                         */
                        
                        double epsilon_pivot = (Q[2][idx_cell_pivot_E] -
                            0.5*Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom]/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                        
                        double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getPressure(
                                &Q[0][idx_cell_pivot_rho],
                                &epsilon_pivot,
                                thermo_properties_ptr);
                        
                        double p = p_pivot;
                        
                        double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getTemperature(
                                &Q[0][idx_cell_pivot_rho],
                                &p_pivot,
                                thermo_properties_ptr);
                        
                        double T = -T_pivot + 2.0*d_bdry_node_isothermal_no_slip_T[node_loc];
                        
                        double rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getDensity(
                                &p,
                                &T,
                                thermo_properties_ptr);
                        
                        double u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                            2.0*d_bdry_node_isothermal_no_slip_vel[node_loc];
                        
                        Q[0][idx_cell_rho] = rho;
                        Q[1][idx_cell_mom] = rho*u;
                        
                        double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                            getInternalEnergyFromTemperature(
                                &rho,
                                &T,
                                thermo_properties_ptr);
                        
                        double E = rho*epsilon + 0.5*Q[1][idx_cell_mom]*Q[1][idx_cell_mom]/rho;
                        
                        Q[2][idx_cell_E] = E;
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
            }
        }
    }
}


/*
 * Function to fill 2d edge boundary values for a patch.
 * Edge locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::fill2dEdgeBoundaryData(
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
    const hier::Patch& patch,
    std::vector<int>& bdry_edge_locs,
    const std::vector<int>& bdry_edge_conds,
    const std::vector<std::vector<double> >& bdry_edge_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(static_cast<int>(conservative_var_data.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT(conservative_var_data[vi]);
    }
    TBOX_ASSERT(static_cast<int>(bdry_edge_locs.size()) <= NUM_2D_EDGES);
    TBOX_ASSERT(*min_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) < NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_conds.size()) == NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_values.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(bdry_edge_values.size()); vi++)
    {
        TBOX_ASSERT(static_cast<int>(bdry_edge_values[vi].size()) ==
                    NUM_2D_EDGES*(conservative_var_data[vi]->getDepth()));
    }
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(2));
    
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY3(*conservative_var_data[vi], patch, ghost_width_to_fill);
    }
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    hier::IntVector num_ghosts(conservative_var_data[0]->getGhostCellWidth());
    for (int vi = 1; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        num_ghosts = hier::IntVector::min(
            num_ghosts,
            conservative_var_data[vi]->getGhostCellWidth());
    }
    
    /*
     * Determine the ghost cell width to fill.
     */
    
    hier::IntVector gcw_to_fill(tbox::Dimension(2));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(2)))
    {
        gcw_to_fill = num_ghosts;
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
    
    const std::vector<hier::BoundaryBox>& edge_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::EDGE2D);
    
    for (int ei = 0; ei < static_cast<int>(edge_bdry.size()); ei++)
    {
        TBOX_ASSERT(edge_bdry[ei].getBoundaryType() == BDRY::EDGE2D);
        
        int edge_loc = edge_bdry[ei].getLocationIndex();
        
        if (std::find(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc) !=
            bdry_edge_locs.end())
        {
            hier::Box fill_box(patch_geom->getBoundaryFillBox(
                edge_bdry[ei],
                interior_box,
                gcw_to_fill));
            
            hier::Index fill_box_lo_idx(fill_box.lower());
            hier::Index fill_box_hi_idx(fill_box.upper());
            
            /*
             * Offset the indices.
             */
            
            fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
            fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
            
            if ((bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP) ||
                (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP))
            {
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<double*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const double*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho = idx_cell_rho;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                        num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                        num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                        num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                        num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                        num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                        num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0];
                            }
                            
                            /*
                             * Set the values for density and momentum.
                             */
                            
                            Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                            Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                2.0*Q[0][idx_cell_pivot_rho]*d_bdry_edge_adiabatic_no_slip_vel[edge_loc*2];
                            Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                2.0*Q[0][idx_cell_pivot_rho]*d_bdry_edge_adiabatic_no_slip_vel[edge_loc*2 + 1];
                            
                            /*
                             * Set the values for total internal energy.
                             */
                            
                            double epsilon_pivot = (Q[3][idx_cell_pivot_E] -
                                0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                     Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom])/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                            
                            double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &Q[0][idx_cell_pivot_rho],
                                    &epsilon_pivot,
                                    thermo_properties_ptr);
                            
                            double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getTemperature(
                                    &Q[0][idx_cell_pivot_rho],
                                    &p_pivot,
                                    thermo_properties_ptr);
                            
                            double T = T_pivot;
                            
                            double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergyFromTemperature(
                                    &Q[0][idx_cell_rho],
                                    &T,
                                    thermo_properties_ptr);
                            
                            double E = Q[0][idx_cell_rho]*epsilon +
                                0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
                else if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho = idx_cell_rho;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                        num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                        num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                        num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                        num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                        num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                        num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0];
                            }
                            
                            /*
                             * Set the values for density, momentum and total internal energy.
                             */
                            
                            double epsilon_pivot = (Q[3][idx_cell_pivot_E] -
                                0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                     Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom])/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                            
                            double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &Q[0][idx_cell_pivot_rho],
                                    &epsilon_pivot,
                                    thermo_properties_ptr);
                            
                            double p = p_pivot;
                            
                            double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getTemperature(
                                    &Q[0][idx_cell_pivot_rho],
                                    &p_pivot,
                                    thermo_properties_ptr);
                            
                            double T = -T_pivot + 2.0*d_bdry_edge_isothermal_no_slip_T[edge_loc];
                            
                            double rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getDensity(
                                    &p,
                                    &T,
                                    thermo_properties_ptr);
                            
                            double u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                2.0*d_bdry_edge_isothermal_no_slip_vel[edge_loc*2];
                            double v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                2.0*d_bdry_edge_isothermal_no_slip_vel[edge_loc*2 + 1];
                            
                            Q[0][idx_cell_rho] = rho;
                            Q[1][idx_cell_mom] = rho*u;
                            Q[2][idx_cell_mom] = rho*v;
                            
                            double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergyFromTemperature(
                                    &rho,
                                    &T,
                                    thermo_properties_ptr);
                            
                            double E = rho*epsilon +
                                0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                    rho;
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
            }
        }
    }
}


/*
 * Function to fill 2d node boundary values for a patch.
 * Node locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::fill2dNodeBoundaryData(
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
    const hier::Patch& patch,
    std::vector<int>& bdry_node_locs,
    const std::vector<int>& bdry_node_conds,
    const std::vector<std::vector<double> >& bdry_edge_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(static_cast<int>(conservative_var_data.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT(conservative_var_data[vi]);
    }
    TBOX_ASSERT(static_cast<int>(bdry_node_locs.size()) <= NUM_2D_NODES);
    TBOX_ASSERT(*min_element(bdry_node_locs.begin(), bdry_node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_node_locs.begin(), bdry_node_locs.end()) < NUM_2D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_2D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_values.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(bdry_edge_values.size()); vi++)
    {
        TBOX_ASSERT(static_cast<int>(bdry_edge_values[vi].size()) ==
                    NUM_2D_EDGES*(conservative_var_data[vi]->getDepth()));
    }
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(2));
    
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY3(*conservative_var_data[vi], patch, ghost_width_to_fill);
    }
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    hier::IntVector num_ghosts(conservative_var_data[0]->getGhostCellWidth());
    for (int vi = 1; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        num_ghosts = hier::IntVector::min(
            num_ghosts,
            conservative_var_data[vi]->getGhostCellWidth());
    }
    
    /*
     * Determine the ghost cell width to fill.
     */
    
    hier::IntVector gcw_to_fill(tbox::Dimension(2));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(2)))
    {
        gcw_to_fill = num_ghosts;
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
    
    const std::vector<hier::BoundaryBox>& node_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::NODE2D);
    
    for (int ni = 0; ni < static_cast<int>(node_bdry.size()); ni++)
    {
        TBOX_ASSERT(node_bdry[ni].getBoundaryType() == BDRY::NODE2D);
        
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
            
            int edge_loc_0 = -1;
            int edge_loc_1 = -1;
            
            switch (node_loc)
            {
                case NODE_BDRY_LOC_2D::XLO_YLO:
                {
                    edge_loc_0 = BDRY_LOC::XLO;
                    edge_loc_1 = BDRY_LOC::YLO;
                    
                    break;
                }
                case NODE_BDRY_LOC_2D::XHI_YLO:
                {
                    edge_loc_0 = BDRY_LOC::XHI;
                    edge_loc_1 = BDRY_LOC::YLO;
                    
                    break;
                }
                case NODE_BDRY_LOC_2D::XLO_YHI:
                {
                    edge_loc_0 = BDRY_LOC::XLO;
                    edge_loc_1 = BDRY_LOC::YHI;
                    
                    break;
                }
                case NODE_BDRY_LOC_2D::XHI_YHI:
                {
                    edge_loc_0 = BDRY_LOC::XHI;
                    edge_loc_1 = BDRY_LOC::YHI;
                    
                    break;
                }
            }
            
            if ((bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP))
            {
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<double*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const double*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho = idx_cell_rho;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc_0 == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc_0 == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            
                            /*
                             * Set the values for density and momentum.
                             */
                            
                            Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                            Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                2.0*Q[0][idx_cell_pivot_rho]*d_bdry_edge_adiabatic_no_slip_vel[edge_loc_0*2];
                            Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                2.0*Q[0][idx_cell_pivot_rho]*d_bdry_edge_adiabatic_no_slip_vel[edge_loc_0*2 + 1];
                            
                            /*
                             * Set the values for total internal energy.
                             */
                            
                            double epsilon_pivot = (Q[3][idx_cell_pivot_E] -
                                0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                     Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom])/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                            
                            double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &Q[0][idx_cell_pivot_rho],
                                    &epsilon_pivot,
                                    thermo_properties_ptr);
                            
                            double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getTemperature(
                                    &Q[0][idx_cell_pivot_rho],
                                    &p_pivot,
                                    thermo_properties_ptr);
                            
                            double T = T_pivot;
                            
                            double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergyFromTemperature(
                                    &Q[0][idx_cell_rho],
                                    &T,
                                    thermo_properties_ptr);
                            
                            double E = Q[0][idx_cell_rho]*epsilon +
                                0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho = idx_cell_rho;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc_1 == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc_1 == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            
                            /*
                             * Set the values for density and momentum.
                             */
                            
                            Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                            Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                2.0*Q[0][idx_cell_pivot_rho]*d_bdry_edge_adiabatic_no_slip_vel[edge_loc_1*2];
                            Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                2.0*Q[0][idx_cell_pivot_rho]*d_bdry_edge_adiabatic_no_slip_vel[edge_loc_1*2 + 1];
                            
                            /*
                             * Set the values for total internal energy.
                             */
                            
                            double epsilon_pivot = (Q[3][idx_cell_pivot_E] -
                                0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                     Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom])/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                            
                            double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &Q[0][idx_cell_pivot_rho],
                                    &epsilon_pivot,
                                    thermo_properties_ptr);
                            
                            double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getTemperature(
                                    &Q[0][idx_cell_pivot_rho],
                                    &p_pivot,
                                    thermo_properties_ptr);
                            
                            double T = T_pivot;
                            
                            double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergyFromTemperature(
                                    &Q[0][idx_cell_rho],
                                    &T,
                                    thermo_properties_ptr);
                            
                            double E = Q[0][idx_cell_rho]*epsilon +
                                0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                Q[0][idx_cell_rho];
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho = idx_cell_rho;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc_0 == BDRY_LOC::XLO)
                            {
                                idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc_0 == BDRY_LOC::XHI)
                            {
                                idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            
                            /*
                             * Set the values for density, momentum and total internal energy.
                             */
                            
                            double epsilon_pivot = (Q[3][idx_cell_pivot_E] -
                                0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                     Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom])/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                            
                            double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &Q[0][idx_cell_pivot_rho],
                                    &epsilon_pivot,
                                    thermo_properties_ptr);
                            
                            double p = p_pivot;
                            
                            double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getTemperature(
                                    &Q[0][idx_cell_pivot_rho],
                                    &p_pivot,
                                    thermo_properties_ptr);
                            
                            double T = -T_pivot + 2.0*d_bdry_edge_isothermal_no_slip_T[edge_loc_0];
                            
                            double rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getDensity(
                                    &p,
                                    &T,
                                    thermo_properties_ptr);
                            
                            double u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                2.0*d_bdry_edge_isothermal_no_slip_vel[edge_loc_0*2];
                            double v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                2.0*d_bdry_edge_isothermal_no_slip_vel[edge_loc_0*2 + 1];
                            
                            Q[0][idx_cell_rho] = rho;
                            Q[1][idx_cell_mom] = rho*u;
                            Q[2][idx_cell_mom] = rho*v;
                            
                            double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergyFromTemperature(
                                    &rho,
                                    &T,
                                    thermo_properties_ptr);
                            
                            double E = rho*epsilon +
                                0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                    rho;
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                (j + num_subghosts_conservative_var[0][1])*
                                    subghostcell_dims_conservative_var[0][0];
                            
                            const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                (j + num_subghosts_conservative_var[1][1])*
                                    subghostcell_dims_conservative_var[1][0];
                            
                            const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                (j + num_subghosts_conservative_var[2][1])*
                                    subghostcell_dims_conservative_var[2][0];
                            
                            int idx_cell_pivot_rho = idx_cell_rho;
                            int idx_cell_pivot_mom = idx_cell_mom;
                            int idx_cell_pivot_E = idx_cell_E;
                            
                            if (edge_loc_1 == BDRY_LOC::YLO)
                            {
                                idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            else if (edge_loc_1 == BDRY_LOC::YHI)
                            {
                                idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0];
                                
                                idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0];
                                
                                idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                        num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0];
                            }
                            
                            /*
                             * Set the values for density, momentum and total internal energy.
                             */
                            
                            double epsilon_pivot = (Q[3][idx_cell_pivot_E] -
                                0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                     Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom])/
                                Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                            
                            double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getPressure(
                                    &Q[0][idx_cell_pivot_rho],
                                    &epsilon_pivot,
                                    thermo_properties_ptr);
                            
                            double p = p_pivot;
                            
                            double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getTemperature(
                                    &Q[0][idx_cell_pivot_rho],
                                    &p_pivot,
                                    thermo_properties_ptr);
                            
                            double T = -T_pivot + 2.0*d_bdry_edge_isothermal_no_slip_T[edge_loc_1];
                            
                            double rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getDensity(
                                    &p,
                                    &T,
                                    thermo_properties_ptr);
                            
                            double u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                2.0*d_bdry_edge_isothermal_no_slip_vel[edge_loc_1*2];
                            double v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                2.0*d_bdry_edge_isothermal_no_slip_vel[edge_loc_1*2 + 1];
                            
                            Q[0][idx_cell_rho] = rho;
                            Q[1][idx_cell_mom] = rho*u;
                            Q[2][idx_cell_mom] = rho*v;
                            
                            double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                getInternalEnergyFromTemperature(
                                    &rho,
                                    &T,
                                    thermo_properties_ptr);
                            
                            double E = rho*epsilon +
                                0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom])/
                                    rho;
                            
                            Q[3][idx_cell_E] = E;
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
            }
        }
    }
}


/*
 * Function to fill 3d face boundary values for a patch.
 * Face locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::fill3dFaceBoundaryData(
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
    const hier::Patch& patch,
    std::vector<int>& bdry_face_locs,
    const std::vector<int>& bdry_face_conds,
    const std::vector<std::vector<double> >& bdry_face_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(static_cast<int>(conservative_var_data.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT(conservative_var_data[vi]);
    }
    TBOX_ASSERT(static_cast<int>(bdry_face_locs.size()) <= NUM_3D_FACES);
    TBOX_ASSERT(*min_element(bdry_face_locs.begin(), bdry_face_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_face_locs.begin(), bdry_face_locs.end()) < NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(bdry_face_conds.size()) == NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(bdry_face_values.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(bdry_face_values.size()); vi++)
    {
        TBOX_ASSERT(static_cast<int>(bdry_face_values[vi].size()) ==
                    NUM_3D_FACES*(conservative_var_data[vi]->getDepth()));
    }
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(3));
    
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY3(*conservative_var_data[vi], patch, ghost_width_to_fill);
    }
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    hier::IntVector num_ghosts(conservative_var_data[0]->getGhostCellWidth());
    for (int vi = 1; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        num_ghosts = hier::IntVector::min(
            num_ghosts,
            conservative_var_data[vi]->getGhostCellWidth());
    }
    
    /*
     * Determine the ghost cell width to fill.
     */
    
    hier::IntVector gcw_to_fill(tbox::Dimension(3));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(3)))
    {
        gcw_to_fill = num_ghosts;
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
    
    const std::vector<hier::BoundaryBox>& face_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::FACE3D);
    
    for (int fi = 0; fi < static_cast<int>(face_bdry.size()); fi++)
    {
        TBOX_ASSERT(face_bdry[fi].getBoundaryType() == BDRY::FACE3D);
        
        int face_loc = face_bdry[fi].getLocationIndex();
        
        if (std::find(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc) !=
            bdry_face_locs.end())
        {
            hier::Box fill_box(patch_geom->getBoundaryFillBox(
                face_bdry[fi],
                interior_box,
                gcw_to_fill));
            
            hier::Index fill_box_lo_idx(fill_box.lower());
            hier::Index fill_box_hi_idx(fill_box.upper());
            
            /*
             * Offset the indices.
             */
            
            fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
            fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
            
            if ((bdry_face_conds[face_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP) ||
                (bdry_face_conds[face_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP))
            {
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<double*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const double*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (bdry_face_conds[face_loc] == BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::XHI)
                                {
                                    idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::YHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::ZHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = T_pivot;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = Q[0][idx_cell_rho]*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove face locations that have boundary conditions identified.
                    bdry_face_locs.erase(std::remove(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc),
                        bdry_face_locs.end());
                }
                else if (bdry_face_conds[face_loc] == BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::XHI)
                                {
                                    idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::YHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc == BDRY_LOC::ZHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double p = p_pivot;
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = -T_pivot + 2.0*d_bdry_face_isothermal_no_slip_T[face_loc];
                                
                                double rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_face_isothermal_no_slip_vel[face_loc*3];
                                double v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_face_isothermal_no_slip_vel[face_loc*3 + 1];
                                double w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_face_isothermal_no_slip_vel[face_loc*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = rho*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove face locations that have boundary conditions identified.
                    bdry_face_locs.erase(std::remove(bdry_face_locs.begin(), bdry_face_locs.end(), face_loc),
                        bdry_face_locs.end());
                }
            }
        }
    }
}


/*
 * Function to fill 3d edge boundary values for a patch.
 * Edge locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::fill3dEdgeBoundaryData(
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
    const hier::Patch& patch,
    std::vector<int>& bdry_edge_locs,
    const std::vector<int>& bdry_edge_conds,
    const std::vector<std::vector<double> >& bdry_face_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(static_cast<int>(conservative_var_data.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT(conservative_var_data[vi]);
    }
    TBOX_ASSERT(static_cast<int>(bdry_edge_locs.size()) <= NUM_3D_EDGES);
    TBOX_ASSERT(*min_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_edge_locs.begin(), bdry_edge_locs.end()) < NUM_3D_EDGES);
    TBOX_ASSERT(static_cast<int>(bdry_edge_conds.size()) == NUM_3D_EDGES);
    TBOX_ASSERT(static_cast<int>(bdry_face_values.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(bdry_face_values.size()); vi++)
    {
        TBOX_ASSERT(static_cast<int>(bdry_face_values[vi].size()) ==
                    NUM_3D_FACES*(conservative_var_data[vi]->getDepth()));
    }
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(3));
    
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY3(*conservative_var_data[vi], patch, ghost_width_to_fill);
    }
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    hier::IntVector num_ghosts(conservative_var_data[0]->getGhostCellWidth());
    for (int vi = 1; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        num_ghosts = hier::IntVector::min(
            num_ghosts,
            conservative_var_data[vi]->getGhostCellWidth());
    }
    
    /*
     * Determine the ghost cell width to fill.
     */
    
    hier::IntVector gcw_to_fill(tbox::Dimension(3));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(3)))
    {
        gcw_to_fill = num_ghosts;
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
    
    const std::vector<hier::BoundaryBox>& edge_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::EDGE3D);
    
    for (int ei = 0; ei < static_cast<int>(edge_bdry.size()); ei++)
    {
        TBOX_ASSERT(edge_bdry[ei].getBoundaryType() == BDRY::EDGE3D);
        
        int edge_loc(edge_bdry[ei].getLocationIndex());
        
        if (std::find(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc) !=
            bdry_edge_locs.end())
        {
            hier::Box fill_box(patch_geom->getBoundaryFillBox(
                edge_bdry[ei],
                interior_box,
                gcw_to_fill));
            
            hier::Index fill_box_lo_idx(fill_box.lower());
            hier::Index fill_box_hi_idx(fill_box.upper());
            
            /*
             * Offset the indices.
             */
            
            fill_box_lo_idx = fill_box_lo_idx - interior_box.lower();
            fill_box_hi_idx = fill_box_hi_idx - interior_box.lower();
            
            int face_loc_0 = -1;
            int face_loc_1 = -1;
            int face_loc_2 = -1;
            
            switch (edge_loc)
            {
                case EDGE_BDRY_LOC_3D::XLO_YLO:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_1 = BDRY_LOC::YLO;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_YLO:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_1 = BDRY_LOC::YLO;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_YHI:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_1 = BDRY_LOC::YHI;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_YHI:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_1 = BDRY_LOC::YHI;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_ZLO:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_ZLO:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_ZHI:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_ZHI:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::YLO_ZLO:
                {
                    face_loc_1 = BDRY_LOC::YLO;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::YHI_ZLO:
                {
                    face_loc_1 = BDRY_LOC::YHI;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::YLO_ZHI:
                {
                    face_loc_1 = BDRY_LOC::YLO;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
                case EDGE_BDRY_LOC_3D::YHI_ZHI:
                {
                    face_loc_1 = BDRY_LOC::YHI;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
            }
            
            if ((bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP) ||
                (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP) ||
                (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP) ||
                (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP) ||
                (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP) ||
                (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP))
            {
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<double*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const double*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_0 == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_0 == BDRY_LOC::XHI)
                                {
                                    idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = T_pivot;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = Q[0][idx_cell_rho]*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
                else if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_1 == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_1 == BDRY_LOC::YHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = T_pivot;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = Q[0][idx_cell_rho]*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
                else if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_2 == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_2 == BDRY_LOC::ZHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = T_pivot;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = Q[0][idx_cell_rho]*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
                else if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_0 == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_0 == BDRY_LOC::XHI)
                                {
                                    idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double p = p_pivot;
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = -T_pivot + 2.0*d_bdry_face_isothermal_no_slip_T[face_loc_0];
                                
                                double rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3];
                                double v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3 + 1];
                                double w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = rho*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
                else if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_1 == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_1 == BDRY_LOC::YHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double p = p_pivot;
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = -T_pivot + 2.0*d_bdry_face_isothermal_no_slip_T[face_loc_1];
                                
                                double rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3];
                                double v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3 + 1];
                                double w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = rho*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
                else if (bdry_edge_conds[edge_loc] == BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_2 == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_2 == BDRY_LOC::ZHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double p = p_pivot;
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = -T_pivot + 2.0*d_bdry_face_isothermal_no_slip_T[face_loc_2];
                                
                                double rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3];
                                double v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3 + 1];
                                double w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = rho*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove edge locations that have boundary conditions identified.
                    bdry_edge_locs.erase(std::remove(bdry_edge_locs.begin(), bdry_edge_locs.end(), edge_loc),
                        bdry_edge_locs.end());
                }
            }
        }
    }
}


/*
 * Function to fill 3d node boundary values for a patch.
 * Node locations that have boundary conditions identified are removed from the container.
 */
void
FlowModelBoundaryUtilitiesSingleSpecies::fill3dNodeBoundaryData(
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_var_data,
    const hier::Patch& patch,
    std::vector<int>& bdry_node_locs,
    const std::vector<int>& bdry_node_conds,
    const std::vector<std::vector<double> >& bdry_face_values,
    const hier::IntVector& ghost_width_to_fill)
{
    TBOX_ASSERT(static_cast<int>(conservative_var_data.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT(conservative_var_data[vi]);
    }
    TBOX_ASSERT(static_cast<int>(bdry_node_locs.size()) <= NUM_3D_NODES);
    TBOX_ASSERT(*min_element(bdry_node_locs.begin(), bdry_node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(bdry_node_locs.begin(), bdry_node_locs.end()) < NUM_3D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_node_conds.size()) == NUM_3D_NODES);
    TBOX_ASSERT(static_cast<int>(bdry_face_values.size()) == 3);
    for (int vi = 0; vi < static_cast<int>(bdry_face_values.size()); vi++)
    {
        TBOX_ASSERT(static_cast<int>(bdry_face_values[vi].size()) ==
                    NUM_3D_FACES*(conservative_var_data[vi]->getDepth()));
    }
    
    TBOX_DIM_ASSERT(ghost_width_to_fill.getDim() == tbox::Dimension(3));
    
    for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        TBOX_ASSERT_OBJDIM_EQUALITY3(*conservative_var_data[vi], patch, ghost_width_to_fill);
    }
    
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    TBOX_ASSERT(patch_geom);
    
    hier::IntVector num_ghosts(conservative_var_data[0]->getGhostCellWidth());
    for (int vi = 1; vi < static_cast<int>(conservative_var_data.size()); vi++)
    {
        num_ghosts = hier::IntVector::min(
            num_ghosts,
            conservative_var_data[vi]->getGhostCellWidth());
    }
    
    /*
     * Determine the ghost cell width to fill.
     */
    
    hier::IntVector gcw_to_fill(tbox::Dimension(3));
    
    // If the ghost fill width is not used, it is set to the ghost cell width of the data.
    if (ghost_width_to_fill == -hier::IntVector::getOne(tbox::Dimension(3)))
    {
        gcw_to_fill = num_ghosts;
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
    
    const std::vector<hier::BoundaryBox>& node_bdry =
        patch_geom->getCodimensionBoundaries(BDRY::NODE3D);
    
    for (int ni = 0; ni < static_cast<int>(node_bdry.size()); ni++)
    {
        TBOX_ASSERT(node_bdry[ni].getBoundaryType() == BDRY::NODE3D);
        
        int node_loc(node_bdry[ni].getLocationIndex());
        
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
            
            int face_loc_0 = -1;
            int face_loc_1 = -1;
            int face_loc_2 = -1;
            
            switch (node_loc)
            {
                case NODE_BDRY_LOC_3D::XLO_YLO_ZLO:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_1 = BDRY_LOC::YLO;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YLO_ZLO:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_1 = BDRY_LOC::YLO;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XLO_YHI_ZLO:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_1 = BDRY_LOC::YHI;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YHI_ZLO:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_1 = BDRY_LOC::YHI;
                    face_loc_2 = BDRY_LOC::ZLO;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XLO_YLO_ZHI:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_1 = BDRY_LOC::YLO;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YLO_ZHI:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_1 = BDRY_LOC::YLO;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XLO_YHI_ZHI:
                {
                    face_loc_0 = BDRY_LOC::XLO;
                    face_loc_1 = BDRY_LOC::YHI;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YHI_ZHI:
                {
                    face_loc_0 = BDRY_LOC::XHI;
                    face_loc_1 = BDRY_LOC::YHI;
                    face_loc_2 = BDRY_LOC::ZHI;
                    
                    break;
                }
            }
            
            if ((bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP) ||
                (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP))
            {
                /*
                 * Get the pointers to the conservative variables.
                 * Get the numbers of ghost cells and the dimensions of the ghost cell boxes of
                 * the conservative variables.
                 */
                
                std::vector<double*> Q;
                Q.reserve(d_num_eqn);
                
                std::vector<hier::IntVector> num_subghosts_conservative_var;
                std::vector<hier::IntVector> subghostcell_dims_conservative_var;
                
                num_subghosts_conservative_var.reserve(3);
                subghostcell_dims_conservative_var.reserve(3);
                
                int count_eqn = 0;
                
                for (int vi = 0; vi < static_cast<int>(conservative_var_data.size()); vi++)
                {
                    int depth = conservative_var_data[vi]->getDepth();
                    
                    for (int di = 0; di < depth; di++)
                    {
                        // If the last element of the conservative variable vector is not in the
                        // system of equations, ignore it.
                        if (count_eqn >= d_num_eqn)
                            break;
                        
                        Q.push_back(conservative_var_data[vi]->getPointer(di));
                        
                        count_eqn++;
                    }
                    
                    num_subghosts_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostCellWidth());
                    subghostcell_dims_conservative_var.push_back(
                        conservative_var_data[vi]->getGhostBox().numberCells());
                }
                
                // Get the thermodynamic properties of the species.
                std::vector<const double*> thermo_properties_ptr;
                thermo_properties_ptr.reserve(static_cast<int> (d_thermo_properties.size()));
                for (int ti = 0; ti < static_cast<int> (d_thermo_properties.size()); ti++)
                {
                    thermo_properties_ptr.push_back(&d_thermo_properties[ti]);
                }
                
                if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_0 == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_0 == BDRY_LOC::XHI)
                                {
                                    idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_0*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = T_pivot;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = Q[0][idx_cell_rho]*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_1 == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_1 == BDRY_LOC::YHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_1*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = T_pivot;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = Q[0][idx_cell_rho]*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_2 == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_2 == BDRY_LOC::ZHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density and momentum.
                                 */
                                
                                Q[0][idx_cell_rho] = Q[0][idx_cell_pivot_rho];
                                Q[1][idx_cell_mom] = -Q[1][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3];
                                Q[2][idx_cell_mom] = -Q[2][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3 + 1];
                                Q[3][idx_cell_mom] = -Q[3][idx_cell_pivot_mom] +
                                    2.0*Q[0][idx_cell_pivot_rho]*d_bdry_face_adiabatic_no_slip_vel[face_loc_2*3 + 2];
                                
                                /*
                                 * Set the values for total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = T_pivot;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &Q[0][idx_cell_rho],
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = Q[0][idx_cell_rho]*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/
                                    Q[0][idx_cell_rho];
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_0 == BDRY_LOC::XLO)
                                {
                                    idx_cell_pivot_rho = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_lo_idx[0] + (fill_box_hi_idx[0] - i) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_0 == BDRY_LOC::XHI)
                                {
                                    idx_cell_pivot_rho = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (interior_box_hi_idx[0] - (i - fill_box_lo_idx[0]) +
                                            num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double p = p_pivot;
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = -T_pivot + 2.0*d_bdry_face_isothermal_no_slip_T[face_loc_0];
                                
                                double rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3];
                                double v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3 + 1];
                                double w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_0*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = rho*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_1 == BDRY_LOC::YLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_lo_idx[1] + (fill_box_hi_idx[1] - j) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_1 == BDRY_LOC::YHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[0][1])*
                                                subghostcell_dims_conservative_var[0][0] +
                                        (k + num_subghosts_conservative_var[0][2])*
                                            subghostcell_dims_conservative_var[0][0]*
                                                subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[1][1])*
                                                subghostcell_dims_conservative_var[1][0] +
                                        (k + num_subghosts_conservative_var[1][2])*
                                            subghostcell_dims_conservative_var[1][0]*
                                                subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (interior_box_hi_idx[1] - (j - fill_box_lo_idx[1]) +
                                            num_subghosts_conservative_var[2][1])*
                                                subghostcell_dims_conservative_var[2][0] +
                                        (k + num_subghosts_conservative_var[2][2])*
                                            subghostcell_dims_conservative_var[2][0]*
                                                subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double p = p_pivot;
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = -T_pivot + 2.0*d_bdry_face_isothermal_no_slip_T[face_loc_1];
                                
                                double rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3];
                                double v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3 + 1];
                                double w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_1*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = rho*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
                else if (bdry_node_conds[node_loc] == BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP)
                {
                    for (int i = fill_box_lo_idx[0]; i <= fill_box_hi_idx[0]; i++)
                    {
                        for (int j = fill_box_lo_idx[1]; j <= fill_box_hi_idx[1]; j++)
                        {
                            for (int k = fill_box_lo_idx[2]; k <= fill_box_hi_idx[2]; k++)
                            {
                                const int idx_cell_rho = (i + num_subghosts_conservative_var[0][0]) +
                                    (j + num_subghosts_conservative_var[0][1])*
                                        subghostcell_dims_conservative_var[0][0] +
                                    (k + num_subghosts_conservative_var[0][2])*
                                        subghostcell_dims_conservative_var[0][0]*
                                            subghostcell_dims_conservative_var[0][1];
                                
                                const int idx_cell_mom = (i + num_subghosts_conservative_var[1][0]) +
                                    (j + num_subghosts_conservative_var[1][1])*
                                        subghostcell_dims_conservative_var[1][0] +
                                    (k + num_subghosts_conservative_var[1][2])*
                                        subghostcell_dims_conservative_var[1][0]*
                                            subghostcell_dims_conservative_var[1][1];
                                
                                const int idx_cell_E = (i + num_subghosts_conservative_var[2][0]) +
                                    (j + num_subghosts_conservative_var[2][1])*
                                        subghostcell_dims_conservative_var[2][0] +
                                    (k + num_subghosts_conservative_var[2][2])*
                                        subghostcell_dims_conservative_var[2][0]*
                                            subghostcell_dims_conservative_var[2][1];
                                
                                int idx_cell_pivot_rho = idx_cell_rho;
                                int idx_cell_pivot_mom = idx_cell_mom;
                                int idx_cell_pivot_E = idx_cell_E;
                                
                                if (face_loc_2 == BDRY_LOC::ZLO)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_lo_idx[2] + (fill_box_hi_idx[2] - k) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                else if (face_loc_2 == BDRY_LOC::ZHI)
                                {
                                    idx_cell_pivot_rho = (i + num_subghosts_conservative_var[0][0]) +
                                        (j + num_subghosts_conservative_var[0][1])*
                                            subghostcell_dims_conservative_var[0][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[0][2])*
                                                subghostcell_dims_conservative_var[0][0]*
                                                    subghostcell_dims_conservative_var[0][1];
                                    
                                    idx_cell_pivot_mom = (i + num_subghosts_conservative_var[1][0]) +
                                        (j + num_subghosts_conservative_var[1][1])*
                                            subghostcell_dims_conservative_var[1][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[1][2])*
                                                subghostcell_dims_conservative_var[1][0]*
                                                    subghostcell_dims_conservative_var[1][1];
                                    
                                    idx_cell_pivot_E = (i + num_subghosts_conservative_var[2][0]) +
                                        (j + num_subghosts_conservative_var[2][1])*
                                            subghostcell_dims_conservative_var[2][0] +
                                        (interior_box_hi_idx[2] - (k - fill_box_lo_idx[2]) +
                                            num_subghosts_conservative_var[2][2])*
                                                subghostcell_dims_conservative_var[2][0]*
                                                    subghostcell_dims_conservative_var[2][1];
                                }
                                
                                /*
                                 * Set the values for density, momentum and total internal energy.
                                 */
                                
                                double epsilon_pivot = (Q[4][idx_cell_pivot_E] -
                                    0.5*(Q[1][idx_cell_pivot_mom]*Q[1][idx_cell_pivot_mom] +
                                         Q[2][idx_cell_pivot_mom]*Q[2][idx_cell_pivot_mom] +
                                         Q[3][idx_cell_pivot_mom]*Q[3][idx_cell_pivot_mom])/
                                    Q[0][idx_cell_pivot_rho])/Q[0][idx_cell_pivot_rho];
                                
                                double p_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getPressure(
                                        &Q[0][idx_cell_pivot_rho],
                                        &epsilon_pivot,
                                        thermo_properties_ptr);
                                
                                double p = p_pivot;
                                
                                double T_pivot = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getTemperature(
                                        &Q[0][idx_cell_pivot_rho],
                                        &p_pivot,
                                        thermo_properties_ptr);
                                
                                double T = -T_pivot + 2.0*d_bdry_face_isothermal_no_slip_T[face_loc_2];
                                
                                double rho = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getDensity(
                                        &p,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double u = -Q[1][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3];
                                double v = -Q[2][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3 + 1];
                                double w = -Q[3][idx_cell_pivot_mom]/Q[0][idx_cell_pivot_rho] +
                                    2.0*d_bdry_edge_isothermal_no_slip_vel[face_loc_2*3 + 2];
                                
                                Q[0][idx_cell_rho] = rho;
                                Q[1][idx_cell_mom] = rho*u;
                                Q[2][idx_cell_mom] = rho*v;
                                Q[3][idx_cell_mom] = rho*w;
                                
                                double epsilon = d_equation_of_state_mixing_rules->getEquationOfState()->
                                    getInternalEnergyFromTemperature(
                                        &rho,
                                        &T,
                                        thermo_properties_ptr);
                                
                                double E = rho*epsilon +
                                    0.5*(Q[1][idx_cell_mom]*Q[1][idx_cell_mom] + Q[2][idx_cell_mom]*Q[2][idx_cell_mom] +
                                         Q[3][idx_cell_mom]*Q[3][idx_cell_mom])/rho;
                                
                                Q[4][idx_cell_E] = E;
                            }
                        }
                    }
                    
                    // Remove node locations that have boundary conditions identified.
                    bdry_node_locs.erase(std::remove(bdry_node_locs.begin(), bdry_node_locs.end(), node_loc),
                        bdry_node_locs.end());
                }
            }
        }
    }
}


void
FlowModelBoundaryUtilitiesSingleSpecies::read1dBdryNodes(
    const boost::shared_ptr<tbox::Database>& input_db,
    std::vector<int>& node_locs,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(1));
    
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
        // Node boundary input required.
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
                boost::shared_ptr<tbox::Database> bdry_loc_db(
                    input_db->getDatabase(bdry_loc_str));
                std::string bdry_cond_str =
                    bdry_loc_db->getString("boundary_condition");
                if (bdry_cond_str == "ADIABATIC_NO_SLIP")
                {
                    node_conds[s] = BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP;
                    
                    readAdiabaticNoSlip(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    node_locs[ni] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ISOTHERMAL_NO_SLIP")
                {
                    node_conds[s] = BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP;
                    
                    readIsothermalNoSlip(
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


void
FlowModelBoundaryUtilitiesSingleSpecies::read2dBdryEdges(
    const boost::shared_ptr<tbox::Database>& input_db,
    std::vector<int>& edge_locs,
    std::vector<int>& edge_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));
    
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(edge_locs.size()) <= NUM_2D_EDGES);
    TBOX_ASSERT(*min_element(edge_locs.begin(), edge_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(edge_locs.begin(), edge_locs.end()) < NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_2D_EDGES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 2; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 2)
    {
        // Edge boundary input required.
        for (int ei = 0; ei < static_cast<int>(edge_locs.size()); ei++)
        {
            int s = edge_locs[ei];
            
            std::string bdry_loc_str;
            switch (s)
            {
                case BDRY_LOC::XLO:
                {
                    bdry_loc_str = "boundary_edge_xlo";
                    break;
                }
                case BDRY_LOC::XHI:
                {
                    bdry_loc_str = "boundary_edge_xhi";
                    break;
                }
                case BDRY_LOC::YLO:
                {
                    bdry_loc_str = "boundary_edge_ylo";
                    break;
                }
                case BDRY_LOC::YHI:
                {
                    bdry_loc_str = "boundary_edge_yhi";
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
                else if (periodic(1) && (s == BDRY_LOC::YLO || s == BDRY_LOC::YHI))
                {
                    need_data_read = false;
                }
            }
            
            if (need_data_read)
            {
                boost::shared_ptr<tbox::Database> bdry_loc_db(
                    input_db->getDatabase(bdry_loc_str));
                std::string bdry_cond_str =
                    bdry_loc_db->getString("boundary_condition");
                if (bdry_cond_str == "ADIABATIC_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP;
                    
                    readAdiabaticNoSlip(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ISOTHERMAL_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP;
                    
                    readIsothermalNoSlip(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
            } // if (need_data_read)
       } // for (int ei = 0 ...
    } // if (num_per_dirs < 2)
    
    // Remove edge locations that have boundary conditions identified.
    edge_locs.erase(std::remove(edge_locs.begin(), edge_locs.end(), BOGUS_BDRY_LOC), edge_locs.end());
}


void
FlowModelBoundaryUtilitiesSingleSpecies::read2dBdryNodes(
    const boost::shared_ptr<tbox::Database>& input_db,
    std::vector<int>& node_locs,
    const std::vector<int>& edge_conds,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(2));
    
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_2D_NODES);
    TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_2D_NODES);
    TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_2D_EDGES);
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_2D_NODES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 2; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 1)
    {
        // Node boundary data required.
        for (int ni = 0; ni < static_cast<int>(node_locs.size()); ni++)
        {
            int s = node_locs[ni];
            
            std::string bdry_loc_str;
            switch (s)
            {
                case NODE_BDRY_LOC_2D::XLO_YLO:
                {
                    bdry_loc_str = "boundary_node_xlo_ylo";
                    break;
                }
                case NODE_BDRY_LOC_2D::XHI_YLO:
                {
                    bdry_loc_str = "boundary_node_xhi_ylo";
                    break;
                }
                case NODE_BDRY_LOC_2D::XLO_YHI:
                {
                    bdry_loc_str = "boundary_node_xlo_yhi";
                    break;
                }
                case NODE_BDRY_LOC_2D::XHI_YHI:
                {
                    bdry_loc_str = "boundary_node_xhi_yhi";
                    break;
                }
                default: NULL_STATEMENT;
            }
            
            boost::shared_ptr<tbox::Database> bdry_loc_db(
                input_db->getDatabase(bdry_loc_str));
            std::string bdry_cond_str =
                bdry_loc_db->getString("boundary_condition");
            if (bdry_cond_str == "XADIABATIC_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YADIABATIC_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "XISOTHERMAL_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YISOTHERMAL_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            
            std::string proper_edge;
            std::string proper_edge_data;
            bool no_edge_data_found = false;
            if ((bdry_cond_str == "XADIABATIC_NO_SLIP") ||
                (bdry_cond_str == "XISOTHERMAL_NO_SLIP"))
            {
                if (s == NODE_BDRY_LOC_2D::XLO_YLO ||
                    s == NODE_BDRY_LOC_2D::XLO_YHI)
                {
                    proper_edge = "XLO";
                    if (bdry_cond_str == "XADIABATIC_NO_SLIP" &&
                        edge_conds[BDRY_LOC::XLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "XISOTHERMAL_NO_SLIP" &&
                        edge_conds[BDRY_LOC::XLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
                else
                {
                    proper_edge = "XHI";
                    if (bdry_cond_str == "XADIABATIC_NO_SLIP" &&
                        edge_conds[BDRY_LOC::XHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "XISOTHERMAL_NO_SLIP" &&
                        edge_conds[BDRY_LOC::XHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
            }
            else if ((bdry_cond_str == "YADIABATIC_NO_SLIP") ||
                     (bdry_cond_str == "YISOTHERMAL_NO_SLIP"))
            {
                if (s == NODE_BDRY_LOC_2D::XLO_YLO ||
                    s == NODE_BDRY_LOC_2D::XHI_YLO)
                {
                    proper_edge = "YLO";
                    if (bdry_cond_str == "YADIABATIC_NO_SLIP" &&
                        edge_conds[BDRY_LOC::YLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "YISOTHREMAL_NO_SLIP" &&
                        edge_conds[BDRY_LOC::YLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
                else
                {
                    proper_edge = "YHI";
                    if (bdry_cond_str == "YADIABATIC_NO_SLIP" &&
                        edge_conds[BDRY_LOC::YHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "YISOTHERMAL_NO_SLIP" &&
                        edge_conds[BDRY_LOC::YHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_edge_data_found = true;
                        proper_edge_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
            }
            
            if (no_edge_data_found)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBoundaryUtilitiesSingleSpecies::read2dBdryNodes()\n"
                    << "Bdry condition '"
                    << bdry_cond_str
                    << "' found for "
                    << bdry_loc_str
                    << "\n but no "
                    << proper_edge_data
                    << " data found for edge "
                    << proper_edge << std::endl);
            }
        } // for (int ni = 0 ...
    } // if (num_per_dirs < 1)
    
    // Remove node locations that have boundary conditions identified.
    node_locs.erase(std::remove(node_locs.begin(), node_locs.end(), BOGUS_BDRY_LOC), node_locs.end());
}


void
FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryFaces(
    const boost::shared_ptr<tbox::Database>& input_db,
    std::vector<int>& face_locs,
    std::vector<int>& face_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(3));
    
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(face_locs.size()) <= NUM_3D_FACES);
    TBOX_ASSERT(*min_element(face_locs.begin(), face_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(face_locs.begin(), face_locs.end()) < NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(face_conds.size()) == NUM_3D_FACES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 3; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 3)
    {
        // Face boundary input required.
        for (int fi = 0; fi < static_cast<int>(face_locs.size()); fi++)
        {
            int s = face_locs[fi];
            
            std::string bdry_loc_str;
            switch (s)
            {
                case BDRY_LOC::XLO:
                {
                    bdry_loc_str = "boundary_face_xlo";
                    break;
                }
                case BDRY_LOC::XHI:
                {
                    bdry_loc_str = "boundary_face_xhi";
                    break;
                }
                case BDRY_LOC::YLO:
                {
                    bdry_loc_str = "boundary_face_ylo";
                    break;
                }
                case BDRY_LOC::YHI:
                {
                    bdry_loc_str = "boundary_face_yhi";
                    break;
                }
                case BDRY_LOC::ZLO:
                {
                    bdry_loc_str = "boundary_face_zlo";
                    break;
                }
                case BDRY_LOC::ZHI:
                {
                    bdry_loc_str = "boundary_face_zhi";
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
                else if (periodic(1) && (s == BDRY_LOC::YLO || s == BDRY_LOC::YHI))
                {
                    need_data_read = false;
                }
                else if (periodic(2) && (s == BDRY_LOC::ZLO || s == BDRY_LOC::ZHI))
                {
                    need_data_read = false;
                }
            }
            
            if (need_data_read)
            {
                boost::shared_ptr<tbox::Database> bdry_loc_db(
                    input_db->getDatabase(bdry_loc_str));
                std::string bdry_cond_str =
                    bdry_loc_db->getString("boundary_condition");
                if (bdry_cond_str == "ADIABATIC_NO_SLIP")
                {
                    face_conds[s] = BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP;
                    
                    readAdiabaticNoSlip(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    face_locs[fi] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ISOTHERMAL_NO_SLIP")
                {
                    face_conds[s] = BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP;
                    
                    readIsothermalNoSlip(
                        bdry_loc_db,
                        bdry_loc_str,
                        s);
                    
                    face_locs[fi] = BOGUS_BDRY_LOC;
                }
            } // if (need_data_read)
        } // for (int fi = 0 ...
    } // if (num_per_dirs < 3)
    
    // Remove face locations that have boundary conditions identified.
    face_locs.erase(std::remove(face_locs.begin(), face_locs.end(), BOGUS_BDRY_LOC), face_locs.end());
}


void
FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryEdges(
    const boost::shared_ptr<tbox::Database>& input_db,
    std::vector<int>& edge_locs,
    const std::vector<int>& face_conds,
    std::vector<int>& edge_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(3));
    
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(edge_locs.size()) <= NUM_3D_EDGES);
    TBOX_ASSERT(*min_element(edge_locs.begin(), edge_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(edge_locs.begin(), edge_locs.end()) < NUM_3D_EDGES);
    TBOX_ASSERT(static_cast<int>(face_conds.size()) == NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(edge_conds.size()) == NUM_3D_EDGES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 3; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 2)
    {
        // edge boundary input required
        for (int ei = 0; ei < static_cast<int>(edge_locs.size()); ei++)
        {
            int s = edge_locs[ei];
            
            std::string bdry_loc_str;
            switch (s)
            {
                case EDGE_BDRY_LOC_3D::YLO_ZLO:
                {
                    bdry_loc_str = "boundary_edge_ylo_zlo";
                    break;
                }
                case EDGE_BDRY_LOC_3D::YHI_ZLO:
                {
                    bdry_loc_str = "boundary_edge_yhi_zlo";
                    break;
                }
                case EDGE_BDRY_LOC_3D::YLO_ZHI:
                {
                    bdry_loc_str = "boundary_edge_ylo_zhi";
                    break;
                }
                case EDGE_BDRY_LOC_3D::YHI_ZHI:
                {
                    bdry_loc_str = "boundary_edge_yhi_zhi";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_ZLO:
                {
                    bdry_loc_str = "boundary_edge_xlo_zlo";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_ZHI:
                {
                    bdry_loc_str = "boundary_edge_xlo_zhi";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_ZLO:
                {
                    bdry_loc_str = "boundary_edge_xhi_zlo";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_ZHI:
                {
                    bdry_loc_str = "boundary_edge_xhi_zhi";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_YLO:
                {
                    bdry_loc_str = "boundary_edge_xlo_ylo";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_YLO:
                {
                    bdry_loc_str = "boundary_edge_xhi_ylo";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XLO_YHI:
                {
                    bdry_loc_str = "boundary_edge_xlo_yhi";
                    break;
                }
                case EDGE_BDRY_LOC_3D::XHI_YHI:
                {
                    bdry_loc_str = "boundary_edge_xhi_yhi";
                    break;
                }
                default: NULL_STATEMENT;
            }
            
            bool need_data_read = false;
            if (num_per_dirs == 0)
            {
                need_data_read = true;
            }
            else if (periodic(0) &&
                     (s == EDGE_BDRY_LOC_3D::YLO_ZLO ||
                      s == EDGE_BDRY_LOC_3D::YHI_ZLO ||
                      s == EDGE_BDRY_LOC_3D::YLO_ZHI ||
                      s == EDGE_BDRY_LOC_3D::YHI_ZHI))
            {
                need_data_read = true;
            }
            else if (periodic(1) &&
                     (s == EDGE_BDRY_LOC_3D::XLO_ZLO ||
                      s == EDGE_BDRY_LOC_3D::XLO_ZHI ||
                      s == EDGE_BDRY_LOC_3D::XHI_ZLO ||
                      s == EDGE_BDRY_LOC_3D::XHI_ZHI))
            {
                need_data_read = true;
            }
            else if (periodic(2) &&
                     (s == EDGE_BDRY_LOC_3D::XLO_YLO ||
                      s == EDGE_BDRY_LOC_3D::XHI_YLO ||
                      s == EDGE_BDRY_LOC_3D::XLO_YHI ||
                      s == EDGE_BDRY_LOC_3D::XHI_YHI))
            {
                need_data_read = true;
            }
            
            if (need_data_read)
            {
                boost::shared_ptr<tbox::Database> bdry_loc_db(
                   input_db->getDatabase(bdry_loc_str));
                
                std::string bdry_cond_str =
                   bdry_loc_db->getString("boundary_condition");
                if (bdry_cond_str == "XADIABATIC_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "YADIABATIC_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ZADIABATIC_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "XISOTHERMAL_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "YISOTHERMAL_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                else if (bdry_cond_str == "ZISOTHERMAL_NO_SLIP")
                {
                    edge_conds[s] = BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP;
                    edge_locs[ei] = BOGUS_BDRY_LOC;
                }
                
                bool ambiguous_type = false;
                if ((bdry_cond_str == "XADIABATIC_NO_SLIP") ||
                    (bdry_cond_str == "XISOTHERMAL_NO_SLIP"))
                {
                    if (s == EDGE_BDRY_LOC_3D::YLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YHI_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YLO_ZHI ||
                        s == EDGE_BDRY_LOC_3D::YHI_ZHI)
                    {
                        ambiguous_type = true;
                    }
                }
                else if ((bdry_cond_str == "YADIABATIC_NO_SLIP") ||
                         (bdry_cond_str == "YISOTHERMAL_NO_SLIP"))
                {
                    if (s == EDGE_BDRY_LOC_3D::XLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::XLO_ZHI ||
                        s == EDGE_BDRY_LOC_3D::XHI_ZLO ||
                        s == EDGE_BDRY_LOC_3D::XHI_ZHI)
                    {
                        ambiguous_type = true;
                    }
                }
                else if ((bdry_cond_str == "ZADIABATIC_NO_SLIP") ||
                         (bdry_cond_str == "ZISOTHERMAL_NO_SLIP"))
                {
                    if (s == EDGE_BDRY_LOC_3D::XLO_YLO ||
                        s == EDGE_BDRY_LOC_3D::XHI_YLO ||
                        s == EDGE_BDRY_LOC_3D::XLO_YHI ||
                        s == EDGE_BDRY_LOC_3D::XHI_YHI)
                    {
                        ambiguous_type = true;
                    }
                }
                if (ambiguous_type)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryEdges()\n"
                        << "Ambiguous bdry condition '"
                        << bdry_cond_str
                        << "' found for "
                        << bdry_loc_str
                        << std::endl);
                }
                
                std::string proper_face;
                std::string proper_face_data;
                bool no_face_data_found = false;
                if ((bdry_cond_str == "XADIABATIC_NO_SLIP") ||
                    (bdry_cond_str == "XISOTHERMAL_NO_SLIP"))
                {
                    if (s == EDGE_BDRY_LOC_3D::XLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::XLO_ZHI ||
                        s == EDGE_BDRY_LOC_3D::XLO_YLO ||
                        s == EDGE_BDRY_LOC_3D::XLO_YHI)
                    {
                        proper_face = "XLO";
                        if (bdry_cond_str == "XADIABATIC_NO_SLIP" &&
                            face_conds[BDRY_LOC::XLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ADIABATIC_NO_SLIP";
                        }
                        if (bdry_cond_str == "XISOTHERMAL_NO_SLIP" &&
                            face_conds[BDRY_LOC::XLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ISOTHERMAL_NO_SLIP";
                        }
                    }
                    else
                    {
                        proper_face = "XHI";
                        if (bdry_cond_str == "XADIABATIC_NO_SLIP" &&
                            face_conds[BDRY_LOC::XHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ADIABATIC_NO_SLIP";
                        }
                        if (bdry_cond_str == "XISOTHERMAL_NO_SLIP" &&
                            face_conds[BDRY_LOC::XHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ISOTHERMAL_NO_SLIP";
                        }
                    }
                }
                else if ((bdry_cond_str == "YADIABATIC_NO_SLIP") ||
                         (bdry_cond_str == "YISOTHERMAL_NO_SLIP"))
                {
                    if (s == EDGE_BDRY_LOC_3D::YLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YLO_ZHI ||
                        s == EDGE_BDRY_LOC_3D::XLO_YLO ||
                        s == EDGE_BDRY_LOC_3D::XHI_YLO)
                    {
                        proper_face = "YLO";
                        if (bdry_cond_str == "YADIABATIC_NO_SLIP" &&
                            face_conds[BDRY_LOC::YLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ADIABATIC_NO_SLIP";
                        }
                        if (bdry_cond_str == "YISOTHERMAL_NO_SLIP" &&
                            face_conds[BDRY_LOC::YLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ISOTHERMAL_NO_SLIP";
                        }
                    }
                    else
                    {
                        proper_face = "YHI";
                        if (bdry_cond_str == "YADIABATIC_NO_SLIP" &&
                            face_conds[BDRY_LOC::YHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ADIABATIC_NO_SLIP";
                        }
                        if (bdry_cond_str == "YISOTHERMAL_NO_SLIP" &&
                            face_conds[BDRY_LOC::YHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ISOTHERMAL_NO_SLIP";
                        }
                    }
                }
                else if ((bdry_cond_str == "ZADIABATIC_NO_SLIP") ||
                         (bdry_cond_str == "ZISOTHERMAL_NO_SLIP"))
                {
                    if (s == EDGE_BDRY_LOC_3D::XLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YHI_ZLO ||
                        s == EDGE_BDRY_LOC_3D::YLO_ZLO ||
                        s == EDGE_BDRY_LOC_3D::XHI_ZLO)
                    {
                        proper_face = "ZLO";
                        if (bdry_cond_str == "ZADIABATIC_NO_SLIP" &&
                            face_conds[BDRY_LOC::ZLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ADIABATIC_NO_SLIP";
                        }
                        if (bdry_cond_str == "ZISOTHERMAL_NO_SLIP" &&
                            face_conds[BDRY_LOC::ZLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ISOTHERMAL_NO_SLIP";
                        }
                    }
                    else
                    {
                        proper_face = "ZHI";
                        if (bdry_cond_str == "ZADIABATIC_NO_SLIP" &&
                            face_conds[BDRY_LOC::ZHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ADIABATIC_NO_SLIP";
                        }
                        if (bdry_cond_str == "ZISOTHERMAL_NO_SLIP" &&
                            face_conds[BDRY_LOC::ZHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                        {
                            no_face_data_found = true;
                            proper_face_data = "ISOTHERMAL_NO_SLIP";
                        }
                    }
                }
                
                if (no_face_data_found)
                {
                    TBOX_ERROR(d_object_name
                        << ": FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryEdges()\n"
                        << "Bdry condition '"
                        << bdry_cond_str
                        << "' found for "
                        << bdry_loc_str
                        << "\n but no "
                        << proper_face_data
                        << " data found for face "
                        << proper_face << std::endl);
                }
            } // if (need_data_read)
        } // for (int ei = 0 ...
    } // if (num_per_dirs < 2)
    
    // Remove edge locations that have boundary conditions identified.
    edge_locs.erase(std::remove(edge_locs.begin(), edge_locs.end(), BOGUS_BDRY_LOC), edge_locs.end());
}


void
FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryNodes(
    const boost::shared_ptr<tbox::Database>& input_db,
    std::vector<int>& node_locs,
    const std::vector<int>& face_conds,
    std::vector<int>& node_conds,
    const hier::IntVector& periodic)
{
    TBOX_DIM_ASSERT(periodic.getDim() == tbox::Dimension(3));
    
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(static_cast<int>(node_locs.size()) <= NUM_3D_NODES);
    TBOX_ASSERT(*min_element(node_locs.begin(), node_locs.end()) >= 0);
    TBOX_ASSERT(*max_element(node_locs.begin(), node_locs.end()) < NUM_3D_NODES);
    TBOX_ASSERT(static_cast<int>(face_conds.size()) == NUM_3D_FACES);
    TBOX_ASSERT(static_cast<int>(node_conds.size()) == NUM_3D_NODES);
    
    int num_per_dirs = 0;
    for (int id = 0; id < 3; id++)
    {
        if (periodic(id))
            ++num_per_dirs;
    }
    
    if (num_per_dirs < 1)
    {
        // Node boundary data required.
        for (int ni = 0; ni < static_cast<int>(node_locs.size()); ni++)
        {
            int s = node_locs[ni];
            
            std::string bdry_loc_str;
            switch (s)
            {
                case NODE_BDRY_LOC_3D::XLO_YLO_ZLO:
                {
                    bdry_loc_str = "boundary_node_xlo_ylo_zlo";
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YLO_ZLO:
                {
                    bdry_loc_str = "boundary_node_xhi_ylo_zlo";
                    break;
                }
                case NODE_BDRY_LOC_3D::XLO_YHI_ZLO:
                {
                    bdry_loc_str = "boundary_node_xlo_yhi_zlo";
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YHI_ZLO:
                {
                    bdry_loc_str = "boundary_node_xhi_yhi_zlo";
                    break;
                }
                case NODE_BDRY_LOC_3D::XLO_YLO_ZHI:
                {
                    bdry_loc_str = "boundary_node_xlo_ylo_zhi";
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YLO_ZHI:
                {
                    bdry_loc_str = "boundary_node_xhi_ylo_zhi";
                    break;
                }
                case NODE_BDRY_LOC_3D::XLO_YHI_ZHI:
                {
                    bdry_loc_str = "boundary_node_xlo_yhi_zhi";
                    break;
                }
                case NODE_BDRY_LOC_3D::XHI_YHI_ZHI:
                {
                    bdry_loc_str = "boundary_node_xhi_yhi_zhi";
                    break;
                }
                default: NULL_STATEMENT;
            }
            boost::shared_ptr<tbox::Database> bdry_loc_db(
                input_db->getDatabase(bdry_loc_str));
            std::string bdry_cond_str =
                bdry_loc_db->getString("boundary_condition");
            if (bdry_cond_str == "XADIABATIC_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::XADIABATIC_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YADIABATIC_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::YADIABATIC_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "ZADIABATIC_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::ZADIABATIC_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "XISOTHERMAL_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::XISOTHERMAL_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "YISOTHERMAL_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::YISOTHERMAL_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            else if (bdry_cond_str == "ZISOTHERMAL_NO_SLIP")
            {
                node_conds[s] = BDRY_COND::FLOW_MODEL::ZISOTHERMAL_NO_SLIP;
                node_locs[ni] = BOGUS_BDRY_LOC;
            }
            
            std::string proper_face;
            std::string proper_face_data;
            bool no_face_data_found = false;
            if ((bdry_cond_str == "XADIABATIC_NO_SLIP") ||
                (bdry_cond_str == "XISOTHERMAL_NO_SLIP"))
            {
                if (s == NODE_BDRY_LOC_3D::XLO_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XLO_YHI_ZLO ||
                    s == NODE_BDRY_LOC_3D::XLO_YLO_ZHI ||
                    s == NODE_BDRY_LOC_3D::XLO_YHI_ZHI)
                {
                    proper_face = "XLO";
                    if (bdry_cond_str == "XADIABATIC_NO_SLIP" &&
                        face_conds[BDRY_LOC::XLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "XISOTHERMAL_NO_SLIP" &&
                        face_conds[BDRY_LOC::XLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
                else
                {
                    proper_face = "XHI";
                    if (bdry_cond_str == "XADIABATIC_NO_SLIP" &&
                        face_conds[BDRY_LOC::XHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "XISOTHERMAL_NO_SLIP" &&
                        face_conds[BDRY_LOC::XHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
            }
            else if ((bdry_cond_str == "YADIABATIC_NO_SLIP") ||
                     (bdry_cond_str == "YISOTHERMAL_NO_SLIP"))
            {
                if (s == NODE_BDRY_LOC_3D::XLO_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XHI_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XLO_YLO_ZHI ||
                    s == NODE_BDRY_LOC_3D::XHI_YLO_ZHI)
                {
                    proper_face = "YLO";
                    if (bdry_cond_str == "YADIABATIC_NO_SLIP" &&
                        face_conds[BDRY_LOC::YLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "YISOTHERMAL_NO_SLIP" &&
                        face_conds[BDRY_LOC::YLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
                else
                {
                    proper_face = "YHI";
                    if (bdry_cond_str == "YADIABATIC_NO_SLIP" &&
                        face_conds[BDRY_LOC::YHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "YISOTHERMAL_NO_SLIP" &&
                        face_conds[BDRY_LOC::YHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
            }
            else if ((bdry_cond_str == "ZADIABATIC_NO_SLIP") ||
                     (bdry_cond_str == "ZISOTHERMAL_NO_SLIP"))
            {
                if (s == NODE_BDRY_LOC_3D::XLO_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XHI_YLO_ZLO ||
                    s == NODE_BDRY_LOC_3D::XLO_YHI_ZLO ||
                    s == NODE_BDRY_LOC_3D::XHI_YHI_ZLO)
                {
                    proper_face = "ZLO";
                    if (bdry_cond_str == "ZADIABATIC_NO_SLIP" &&
                        face_conds[BDRY_LOC::ZLO] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "ZISOTHERMAL_NO_SLIP" &&
                        face_conds[BDRY_LOC::ZLO] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
                else
                {
                    proper_face = "ZHI";
                    if (bdry_cond_str == "ZADIABATIC_NO_SLIP" &&
                        face_conds[BDRY_LOC::ZHI] != BDRY_COND::FLOW_MODEL::ADIABATIC_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ADIABATIC_NO_SLIP";
                    }
                    if (bdry_cond_str == "ZISOTHERMAL_NO_SLIP" &&
                        face_conds[BDRY_LOC::ZHI] != BDRY_COND::FLOW_MODEL::ISOTHERMAL_NO_SLIP)
                    {
                        no_face_data_found = true;
                        proper_face_data = "ISOTHERMAL_NO_SLIP";
                    }
                }
            }
            
            if (no_face_data_found)
            {
                TBOX_ERROR(d_object_name
                    << ": FlowModelBoundaryUtilitiesSingleSpecies::read3dBdryNodes()\n"
                    << "Bdry condition '"
                    << bdry_cond_str
                    << "' found for "
                    << bdry_loc_str
                    << "\n but no "
                    << proper_face_data
                    << " data found for face "
                    << proper_face << std::endl);
            }
        } // for (int ni = 0 ...
    } // if (num_per_dirs < 1)
    
    // Remove node locations that have boundary conditions identified.
    node_locs.erase(std::remove(node_locs.begin(), node_locs.end(), BOGUS_BDRY_LOC), node_locs.end());
}


void
FlowModelBoundaryUtilitiesSingleSpecies::readAdiabaticNoSlip(
    const boost::shared_ptr<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    
    std::vector<double> data_vel;
    
    if (db->keyExists("velocity"))
    {
        data_vel = db->getDoubleVector("velocity");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::readAdiabaticNoSlip()\n"
            << "'velocity' entry missing from '"
            << db_name
            << "' input database."
            << std::endl);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        if (static_cast<int>(data_vel.size()) != 1)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBoundaryUtilitiesSingleSpecies::readAdiabaticNoSlip()\n"
                << "'velocity' entry from '"
                << db_name
                << "' input database has incorrect size."
                << std::endl);
        }
        
        d_bdry_node_adiabatic_no_slip_vel[bdry_location_index] = data_vel[0];
    }
    else if (d_dim == tbox::Dimension(2))
    {
        if (static_cast<int>(data_vel.size()) != 2)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBoundaryUtilitiesSingleSpecies::readAdiabaticNoSlip()\n"
                << "'velocity' entry from '"
                << db_name
                << "' input database has incorrect size."
                << std::endl);
        }
        
        d_bdry_edge_adiabatic_no_slip_vel[bdry_location_index*2] = data_vel[0];
        d_bdry_edge_adiabatic_no_slip_vel[bdry_location_index*2 + 1] = data_vel[1];
    }
    else if (d_dim == tbox::Dimension(3))
    {
        if (static_cast<int>(data_vel.size()) != 3)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBoundaryUtilitiesSingleSpecies::readAdiabaticNoSlip()\n"
                << "'velocity' entry from '"
                << db_name
                << "' input database has incorrect size."
                << std::endl);
        }
        
        d_bdry_face_adiabatic_no_slip_vel[bdry_location_index*3] = data_vel[0];
        d_bdry_face_adiabatic_no_slip_vel[bdry_location_index*3 + 1] = data_vel[1];
        d_bdry_face_adiabatic_no_slip_vel[bdry_location_index*3 + 2] = data_vel[2];
    }
}


void
FlowModelBoundaryUtilitiesSingleSpecies::readIsothermalNoSlip(
    const boost::shared_ptr<tbox::Database>& db,
    std::string& db_name,
    int bdry_location_index)
{
    TBOX_ASSERT(db);
    TBOX_ASSERT(!db_name.empty());
    
    double data_T = 0.0;
    std::vector<double> data_vel;
    
    if (db->keyExists("temperature"))
    {
        data_T = db->getDouble("temperature");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::readIsothermalNoSlip()\n"
            << "'temperature' entry missing from '"
            << db_name
            << "' input database."
            << std::endl);
    }
    
    if (db->keyExists("velocity"))
    {
        data_vel = db->getDoubleVector("velocity");
    }
    else
    {
        TBOX_ERROR(d_object_name
            << ": FlowModelBoundaryUtilitiesSingleSpecies::readIsothermalNoSlip()\n"
            << "'velocity' entry missing from '"
            << db_name
            << "' input database."
            << std::endl);
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        if (static_cast<int>(data_vel.size()) != 1)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBoundaryUtilitiesSingleSpecies::readIsothermalNoSlip()\n"
                << "'velocity' entry from '"
                << db_name
                << "' input database has incorrect size."
                << std::endl);
        }
        
        d_bdry_node_isothermal_no_slip_T[bdry_location_index] = data_T;
        d_bdry_node_isothermal_no_slip_vel[bdry_location_index] = data_vel[0];
    }
    else if (d_dim == tbox::Dimension(2))
    {
        if (static_cast<int>(data_vel.size()) != 2)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBoundaryUtilitiesSingleSpecies::readIsothermalNoSlip()\n"
                << "'velocity' entry from '"
                << db_name
                << "' input database has incorrect size."
                << std::endl);
        }
        
        d_bdry_edge_isothermal_no_slip_T[bdry_location_index] = data_T;
        d_bdry_edge_isothermal_no_slip_vel[bdry_location_index*2] = data_vel[0];
        d_bdry_edge_isothermal_no_slip_vel[bdry_location_index*2 + 1] = data_vel[1];
    }
    else if (d_dim == tbox::Dimension(3))
    {
        if (static_cast<int>(data_vel.size()) != 3)
        {
            TBOX_ERROR(d_object_name
                << ": FlowModelBoundaryUtilitiesSingleSpecies::readIsothermalNoSlip()\n"
                << "'velocity' entry from '"
                << db_name
                << "' input database has incorrect size."
                << std::endl);
        }
        
        d_bdry_face_isothermal_no_slip_T[bdry_location_index] = data_T;
        d_bdry_face_isothermal_no_slip_vel[bdry_location_index*3] = data_vel[0];
        d_bdry_face_isothermal_no_slip_vel[bdry_location_index*3 + 1] = data_vel[1];
        d_bdry_face_isothermal_no_slip_vel[bdry_location_index*3 + 2] = data_vel[2];
    }
}
