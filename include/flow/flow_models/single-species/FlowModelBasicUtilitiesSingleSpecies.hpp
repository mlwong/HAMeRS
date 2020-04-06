#ifndef FLOW_MODEL_BASIC_UTILITIES_SINGLE_SPECIES_HPP
#define FLOW_MODEL_BASIC_UTILITIES_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModelBasicUtilities.hpp"

class FlowModelBasicUtilitiesSingleSpecies: public FlowModelBasicUtilities
{
    public:
        FlowModelBasicUtilitiesSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules):
                FlowModelBasicUtilities(
                    object_name,
                    dim,
                    grid_geometry,
                    num_species),
            d_equation_of_state_mixing_rules(equation_of_state_mixing_rules)
        {}
        
        ~FlowModelBasicUtilitiesSingleSpecies() {}
        
        /*
         * Check whether the given cell conservative variables are within the bounds.
         */
        void
        checkCellDataOfConservativeVariablesBounded(
            boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables);
        
        /*
         * Check whether the given side conservative variables are within the bounds.
         */
        void
        checkSideDataOfConservativeVariablesBounded(
            boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables);
        
        /*
         * Check whether the given cell primitive variables are within the bounds.
         */
        void
        checkCellDataOfPrimitiveVariablesBounded(
            boost::shared_ptr<pdat::CellData<int> >& bounded_flag,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables);
        
        /*
         * Check whether the given side primitive variables are within the bounds.
         */
        void
        checkSideDataOfPrimitiveVariablesBounded(
            boost::shared_ptr<pdat::SideData<int> >& bounded_flag,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables);
        
    private:
        /*
         * boost::shared_ptr to EquationOfStateMixingRules.
         */
        const boost::shared_ptr<EquationOfStateMixingRules>
            d_equation_of_state_mixing_rules;
        
};

#endif /* FLOW_MODEL_BASIC_UTILITIES_SINGLE_SPECIES_HPP */
