#ifndef FLOW_MODEL_BASIC_UTILITIES_SINGLE_SPECIES_HPP
#define FLOW_MODEL_BASIC_UTILITIES_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModelBasicUtilities.hpp"

class FlowModelBasicUtilitiesSingleSpecies: public FlowModelBasicUtilities
{
    public:
        FlowModelBasicUtilitiesSingleSpecies(
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
                    2 + dim.getValue()),
            d_equation_of_state_mixing_rules(equation_of_state_mixing_rules)
        {}
        
        ~FlowModelBasicUtilitiesSingleSpecies() {}
        
        /*
         * Convert conservative variables to primitive variables.
         */
        void
        convertConservativeVariablesToPrimitiveVariables(
            const std::vector<const double*>& conservative_variables,
            const std::vector<double*>& primitive_variables);
        
        /*
         * Convert conservative variables to primitive variables.
         */
        void
        convertConservativeVariablesToPrimitiveVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables);
        
        /*
         * Convert primitive variables to conservative variables.
         */
        void
        convertPrimitiveVariablesToConservativeVariables(
            const std::vector<const double*>& primitive_variables,
            const std::vector<double*>& conservative_variables);
        
        /*
         * Convert primitive variables to conservative variables.
         */
        void
        convertPrimitiveVariablesToConservativeVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables);
        
        /*
         * Check whether the given cell conservative variables are within the bounds.
         */
        void
        checkCellDataOfConservativeVariablesBounded(
            HAMERS_SHARED_PTR<pdat::CellData<int> >& bounded_flag,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables);
        
        /*
         * Check whether the given side conservative variables are within the bounds.
         */
        void
        checkSideDataOfConservativeVariablesBounded(
            HAMERS_SHARED_PTR<pdat::SideData<int> >& bounded_flag,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables);
        
        /*
         * Check whether the given cell primitive variables are within the bounds.
         */
        void
        checkCellDataOfPrimitiveVariablesBounded(
            HAMERS_SHARED_PTR<pdat::CellData<int> >& bounded_flag,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& primitive_variables);
        
        /*
         * Check whether the given side primitive variables are within the bounds.
         */
        void
        checkSideDataOfPrimitiveVariablesBounded(
            HAMERS_SHARED_PTR<pdat::SideData<int> >& bounded_flag,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables);
        
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
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& projection_variables);
        
        /*
         * Compute the side data of the projection variables for transformation between primitive variables and
         * characteristic variables.
         */
        void
        computeSideDataOfProjectionVariablesForPrimitiveVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& projection_variables);
        
        /*
         * Compute the side data of characteristic variables from conservative variables.
         */
        void
        computeSideDataOfCharacteristicVariablesFromConservativeVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& conservative_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& projection_variables,
            const int& idx_offset);
        
        /*
         * Compute the side data of characteristic variables from primitive variables.
         */
        void
        computeSideDataOfCharacteristicVariablesFromPrimitiveVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::CellData<double> > >& primitive_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& projection_variables,
            const int& idx_offset);
        
        /*
         * Compute the side data of conservative variables from characteristic variables.
         */
        void
        computeSideDataOfConservativeVariablesFromCharacteristicVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& conservative_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& projection_variables);
        
        /*
         * Compute the side data of primitive variables from characteristic variables.
         */
        void
        computeSideDataOfPrimitiveVariablesFromCharacteristicVariables(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& primitive_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& characteristic_variables,
            const std::vector<HAMERS_SHARED_PTR<pdat::SideData<double> > >& projection_variables);
        
    private:
        /*
         * HAMERS_SHARED_PTR to EquationOfStateMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfStateMixingRules>
            d_equation_of_state_mixing_rules;
        
};

#endif /* FLOW_MODEL_BASIC_UTILITIES_SINGLE_SPECIES_HPP */
