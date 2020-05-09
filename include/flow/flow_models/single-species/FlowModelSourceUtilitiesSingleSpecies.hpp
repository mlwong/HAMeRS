#ifndef FLOW_MODEL_SOURCE_UTILITIES_SINGLE_SPECIES_HPP
#define FLOW_MODEL_SOURCE_UTILITIES_SINGLE_SPECIES_HPP

#include "flow/flow_models/FlowModelSourceUtilities.hpp"

class FlowModelSourceUtilitiesSingleSpecies: public FlowModelSourceUtilities
{
    public:
        FlowModelSourceUtilitiesSingleSpecies(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<tbox::Database>& flow_model_db,
            const boost::shared_ptr<EquationOfStateMixingRules> equation_of_state_mixing_rules);
        
        ~FlowModelSourceUtilitiesSingleSpecies() {}
        
        /*
         * Register the required variables for the computation of source terms in the registered patch.
         */
        void
        registerDerivedVariablesForSourceTerms(
            const hier::IntVector& num_subghosts);
        
        /*
         * Allocate memory for cell data of different registered derived variables related to this
         * class in the registered patch.
         */
        void allocateMemoryForDerivedCellData();
        
        /*
         * Clear cell data of different derived variables related to this class in the registered patch.
         */
        void clearCellData();
        
        /*
         * Compute cell data of different registered derived variables related to this class.
         */
        void computeDerivedCellData();
        
        /*
         * Compute all source terms on a patch.
         */
        void
        computeSourceTermsOnPatch(
            const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
            const double time,
            const double dt,
            const int RK_step_number);
        
        /*
         * Put the characteristics of this class into the restart database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
    private:
        /*
         * Whether there is gravity.
         */
        bool d_has_gravity;
        
        /*
         * Gravity vector.
         */
        std::vector<double> d_gravity;
        
        /*
         * boost::shared_ptr to EquationOfStateMixingRules.
         */
        const boost::shared_ptr<EquationOfStateMixingRules>
            d_equation_of_state_mixing_rules;
        
};

#endif /* FLOW_MODEL_SOURCE_UTILITIES_SINGLE_SPECIES_HPP */
