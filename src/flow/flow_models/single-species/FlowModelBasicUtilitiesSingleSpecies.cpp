#include "flow/flow_models/single-species/FlowModelBasicUtilitiesSingleSpecies.hpp"

/*
 * Check whether the given cell conservative variables are within the bounds.
 */
void
FlowModelBasicUtilitiesSingleSpecies::checkCellDataOfConservativeVariablesBounded(
    boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables)
{
}


/*
 * Check whether the given side conservative variables are within the bounds.
 */
void
FlowModelBasicUtilitiesSingleSpecies::checkSideDataOfConservativeVariablesBounded(
    boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables)
{
}


/*
 * Check whether the given cell primitive variables are within the bounds.
 */
void
FlowModelBasicUtilitiesSingleSpecies::checkCellDataOfPrimitiveVariablesBounded(
    boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables)
{
}


/*
 * Check whether the given side primitive variables are within the bounds.
 */
void
FlowModelBasicUtilitiesSingleSpecies::checkSideDataOfPrimitiveVariablesBounded(
    boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables)
{
}
