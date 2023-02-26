#ifndef FLOW_MODEL_SOURCE_UTILITIES_FIVE_EQN_ALLAIRE_HPP
#define FLOW_MODEL_SOURCE_UTILITIES_FIVE_EQN_ALLAIRE_HPP

#include "flow/flow_models/FlowModelSourceUtilities.hpp"

class FlowModelSourceUtilitiesFiveEqnAllaire: public FlowModelSourceUtilities
{
    public:
        FlowModelSourceUtilitiesFiveEqnAllaire(
            const std::string& object_name,
            const std::string& project_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db,
            const HAMERS_SHARED_PTR<EquationOfStateMixingRules> equation_of_state_mixing_rules);
        
        ~FlowModelSourceUtilitiesFiveEqnAllaire() {}
        
        /*
         * Register the required variables for the computation of source terms in the registered patch.
         */
        void
        registerDerivedVariablesForSourceTerms(
            const hier::IntVector& num_subghosts);
        
        /*
         * Register the required variables for the computation of local stable time increment for
         * source terms in the registered patch.
         */
        void
        registerDerivedVariablesForSourceTermsStableDt(
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
            const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_source,
            const double time,
            const double dt,
            const int RK_step_number);
        
        /*
         * Get local stable time increment for source terms.
         */
        double
        getStableDtOnPatch();
        
        /*
         * Put the characteristics of this class into the restart database.
         */
        void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
    private:
        /*
         * Set the number of sub-ghost cells of a variable.
         * This function can be called recursively if the variables are computed recursively.
         */
        void
        setNumberOfSubGhosts(
            const hier::IntVector& num_subghosts,
            const std::string& variable_name,
            const std::string& parent_variable_name);
        
        /*
         * Set the ghost boxes of derived cell variables.
         */
        void
        setDerivedCellVariableGhostBoxes();
        
        /*
         * Compute the cell data of Gruneisen parameter in the registered patch.
         */
        void computeCellDataOfGruneisenParameter();
        
        /*
         * Whether there is gravity.
         */
        bool d_has_gravity;
        
        /*
         * Gravity vector.
         */
        std::vector<Real> d_gravity;
        
        /*
         * Number of sub-ghost cells of derived cell data related to this class.
         */
        hier::IntVector d_num_subghosts_gruneisen_parameter;
        
        /*
         * Boxes with sub-ghost cells of derived cell data related to this class.
         */
        hier::Box d_subghost_box_gruneisen_parameter;
        
        /*
         * Dimensions of boxes with sub-ghost cells of derived cell data related to this class.
         */
        hier::IntVector d_subghostcell_dims_gruneisen_parameter;
        
        /*
         * HAMERS_SHARED_PTR to derived cell data related to this class.
         */
        HAMERS_SHARED_PTR<pdat::CellData<Real> > d_data_gruneisen_parameter;
        
        /*
         * Whether derived cell data related to this class is computed.
         */
        bool d_cell_data_computed_gruneisen_parameter;
        
        /*
         * HAMERS_SHARED_PTR to EquationOfStateMixingRules.
         */
        const HAMERS_SHARED_PTR<EquationOfStateMixingRules>
            d_equation_of_state_mixing_rules;
        
};

#endif /* FLOW_MODEL_SOURCE_UTILITIES_FIVE_EQN_ALLAIRE_HPP */
