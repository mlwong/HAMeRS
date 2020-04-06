#include "flow/flow_models/four-eqn_conservative/FlowModelBasicUtilitiesFourEqnConservative.hpp"

/*
 * Check whether the given cell conservative variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::checkCellDataOfConservativeVariablesBounded(
    boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables)
{
}


/*
 * Check whether the given side conservative variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::checkSideDataOfConservativeVariablesBounded(
    boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables)
{
}


/*
 * Check whether the given cell primitive variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::checkCellDataOfPrimitiveVariablesBounded(
    boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables)
{
}


/*
 * Check whether the given side primitive variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFourEqnConservative::checkSideDataOfPrimitiveVariablesBounded(
    boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables)
{
}
