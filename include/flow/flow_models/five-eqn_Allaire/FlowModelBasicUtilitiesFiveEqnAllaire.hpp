#ifndef FLOW_MODEL_BASIC_UTILITIES_FIVE_EQN_ALLAIRE_HPP
#define FLOW_MODEL_BASIC_UTILITIES_FIVE_EQN_ALLAIRE_HPP

#include "flow/flow_models/FlowModelBasicUtilities.hpp"

class FlowModelBasicUtilitiesFiveEqnAllaire: public FlowModelBasicUtilities
{
    public:
        FlowModelBasicUtilitiesFiveEqnAllaire(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_SHARED_PTR<EquationOfStateMixingRules> equation_of_state_mixing_rules):
                FlowModelBasicUtilities(
                    object_name,
                    dim,
                    grid_geometry,
                    num_species,
                    dim.getValue() + 2*num_species),
            d_equation_of_state_mixing_rules(equation_of_state_mixing_rules)
        {
            // Set the bounds for the variables.
            d_Y_bound_lo = Real(-0.001);
            d_Y_bound_up = Real(1.001);
            d_Z_bound_lo = Real(-1000.0);
            d_Z_bound_up = Real(1000.0);
        }
        
        ~FlowModelBasicUtilitiesFiveEqnAllaire() {}
        
        /*
         * Convert conservative variables to primitive variables.
         */
        void
        convertConservativeVariablesToPrimitiveVariables(
            const std::vector<const Real*>& conservative_variables,
            const std::vector<Real*>& primitive_variables);
        
        /*
         * Convert conservative variables to primitive variables.
         */
        void
        convertConservativeVariablesToPrimitiveVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables);
        
        /*
         * Convert primitive variables to conservative variables.
         */
        void
        convertPrimitiveVariablesToConservativeVariables(
            const std::vector<const Real*>& primitive_variables,
            const std::vector<Real*>& conservative_variables);
        
        /*
         * Convert primitive variables to conservative variables.
         */
        void
        convertPrimitiveVariablesToConservativeVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables);
        
        /*
         * Check whether the given cell conservative variables are within the bounds.
         */
        void
        checkCellDataOfConservativeVariablesBounded(
            HAMERS_SHARED_PTR<pdat::CellData<int> >& bounded_flag,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_variables);
        
        /*
         * Check whether the given side conservative variables are within the bounds.
         */
        void
        checkSideDataOfConservativeVariablesBounded(
            HAMERS_SHARED_PTR<pdat::SideData<int> >& bounded_flag,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables);
        
        /*
         * Check whether the given cell primitive variables are within the bounds.
         */
        void
        checkCellDataOfPrimitiveVariablesBounded(
            HAMERS_SHARED_PTR<pdat::CellData<int> >& bounded_flag,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& primitive_variables);
        
        /*
         * Register the required derived variables for transformation between conservative
         * variables and characteristic variables.
         */
        void
        registerDerivedVariablesForCharacteristicProjectionOfConservativeVariables(
            const hier::IntVector& num_subghosts,
            const AVERAGING::TYPE& averaging_type);
        
        /*
         * Register the required derived variables for transformation between primitive variables
         * and characteristic variables.
         */
        void
        registerDerivedVariablesForCharacteristicProjectionOfPrimitiveVariables(
            const hier::IntVector& num_subghosts,
            const AVERAGING::TYPE& averaging_type);
        
        /*
         * Check whether the given side primitive variables are within the bounds.
         */
        void
        checkSideDataOfPrimitiveVariablesBounded(
            HAMERS_SHARED_PTR<pdat::SideData<int> >& bounded_flag,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables);
        
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
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables);
        
        /*
         * Compute the side data of the projection variables for transformation between primitive variables and
         * characteristic variables.
         */
        void
        computeSideDataOfProjectionVariablesForPrimitiveVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables);
        
        /*
         * Compute the side data of characteristic variables from conservative variables.
         */
        void
        computeSideDataOfCharacteristicVariablesFromConservativeVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& characteristic_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& conservative_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables,
            const int& idx_offset);
        
        /*
         * Compute the side data of characteristic variables from primitive variables.
         */
        void
        computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& characteristic_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<Real> > >& primitive_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables,
            const int& idx_offset);
        
        /*
         * Compute the side data of conservative variables from characteristic variables.
         */
        void
        computeSideDataOfConservativeVariablesFromCharacteristicVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& conservative_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& characteristic_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables);
        
        /*
         * Compute the side data of primitive variables from characteristic variables.
         */
        void
        computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& primitive_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& characteristic_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& projection_variables);
        
    private:
        /*
         * Upper and lower bounds on variables.
         */
        Real d_Y_bound_lo;
        Real d_Y_bound_up;
        Real d_Z_bound_lo;
        Real d_Z_bound_up;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfStateMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfStateMixingRules>
            d_equation_of_state_mixing_rules;
        
};

#endif /* FLOW_MODEL_BASIC_UTILITIES_FIVE_EQN_ALLAIRE_HPP */
