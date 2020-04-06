#include "flow/flow_models/five-eqn_Allaire/FlowModelBasicUtilitiesFiveEqnAllaire.hpp"

/*
 * Check whether the given cell conservative variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::checkCellDataOfConservativeVariablesBounded(
    boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables)
{
}


/*
 * Check whether the given side conservative variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::checkSideDataOfConservativeVariablesBounded(
    boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables)
{
}


/*
 * Check whether the given cell primitive variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::checkCellDataOfPrimitiveVariablesBounded(
    boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables)
{
}


/*
 * Check whether the given side primitive variables are within the bounds.
 */
void
FlowModelBasicUtilitiesFiveEqnAllaire::checkSideDataOfPrimitiveVariablesBounded(
    boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
    const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables)
{
}
