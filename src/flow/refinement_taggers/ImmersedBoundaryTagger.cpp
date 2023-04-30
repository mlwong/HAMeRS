#include "flow/refinement_taggers/ImmersedBoundaryTagger.hpp"

#include "util/immersed_boundaries/ImmersedBoundaries.hpp"

ImmersedBoundaryTagger::ImmersedBoundaryTagger(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
    const hier::IntVector& num_cells_buffer):
        d_object_name(object_name),
        d_dim(dim),
        d_grid_geometry(grid_geometry),
        d_num_cells_buffer(num_cells_buffer),
        d_num_immersed_boundary_tagger_ghosts(d_num_cells_buffer)
{
}


/*
 * Print all characteristics of the immersed boundary tagger class.
 */
void
ImmersedBoundaryTagger::printClassData(std::ostream& os) const
{
    os << "\nPrint GradientTagger object..."
       << std::endl;
    
    os << std::endl;
    
    os << "d_num_cells_buffer = " << d_num_cells_buffer;
}


// /*
//  * Put the characteristics of the immersed boundary tagger class into the restart database.
//  */
// void
// ImmersedBoundaryTagger::putToRestart(
//     const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const
// {
//     NULL_USE(restart_db);
// }


/*
 * Tag cells on a patch for refinement using gradient sensors.
 */
void
ImmersedBoundaryTagger::tagCellsOnPatch(
    hier::Patch& patch,
    const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags,
    const HAMERS_SHARED_PTR<hier::VariableContext>& data_context)
{
    
}
