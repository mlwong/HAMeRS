#ifndef FLOW_MODEL_BASIC_UTILITIES_FIVE_EQN_ALLAIRE_HPP
#define FLOW_MODEL_BASIC_UTILITIES_FIVE_EQN_ALLAIRE_HPP

#include "flow/flow_models/FlowModelBasicUtilities.hpp"

class FlowModelBasicUtilitiesFiveEqnAllaire: public FlowModelBasicUtilities
{
    public:
        FlowModelBasicUtilitiesFiveEqnAllaire(
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
        
        ~FlowModelBasicUtilitiesFiveEqnAllaire() {}
        
        /*
         * Convert conservative variables to primitive variables.
         */
        void
        convertConservativeVariablesToPrimitiveVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables);
        
        /*
         * Convert primitive variables to conservative variables.
         */
        void
        convertPrimitiveVariablesToConservativeVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables);
        
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
        
        /*
         * Get the number of projection variables for transformation between conservative
         * variables and characteristic variables.
         */
        int
        getNumberOfProjectionVariablesForConservativeVariables() const;
        
        /*
         * Get the number of projection variables for transformation between primitive variables
         * and characteristic variables.
         */
        int
        getNumberOfProjectionVariablesForPrimitiveVariables() const;
        
        /*
         * Compute the side data of the projection variables for transformation between conservative variables and
         * characteristic variables.
         */
        void
        computeSideDataOfProjectionVariablesForConservativeVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables);
        
        /*
         * Compute the side data of the projection variables for transformation between primitive variables and
         * characteristic variables.
         */
        void
        computeSideDataOfProjectionVariablesForPrimitiveVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables);
        
        /*
         * Compute the side data of characteristic variables from conservative variables.
         */
        void
        computeSideDataOfCharacteristicVariablesFromConservativeVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables,
            const int& idx_offset);
        
        /*
         * Compute the side data of characteristic variables from primitive variables.
         */
        void
        computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& primitive_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables,
            const int& idx_offset);
        
        /*
         * Compute the side data of conservative variables from characteristic variables.
         */
        void
        computeSideDataOfConservativeVariablesFromCharacteristicVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& conservative_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables);
        
        /*
         * Compute the side data of primitive variables from characteristic variables.
         */
        void
        computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            std::vector<boost::shared_ptr<pdat::SideData<double> > >& primitive_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<boost::shared_ptr<pdat::SideData<double> > >& projection_variables);
        
    private:
        /*
         * boost::shared_ptr to EquationOfStateMixingRules.
         */
        const boost::shared_ptr<EquationOfStateMixingRules>
            d_equation_of_state_mixing_rules;
        
};

#endif /* FLOW_MODEL_BASIC_UTILITIES_FIVE_EQN_ALLAIRE_HPP */
