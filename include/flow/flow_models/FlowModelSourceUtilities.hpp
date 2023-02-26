#ifndef FLOW_MODEL_SOURCE_UTILITIES_HPP
#define FLOW_MODEL_SOURCE_UTILITIES_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModel.hpp"
#include "flow/flow_models/FlowModelSpecialSourceTerms.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellData.h"

#include <string>

class FlowModel;

class FlowModelSourceUtilities
{
    public:
        FlowModelSourceUtilities(
            const std::string& object_name,
            const std::string& project_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<tbox::Database>& flow_model_db);
        
        virtual ~FlowModelSourceUtilities() {}
        
        /*
         * Check whether there are any source terms.
         */
        bool hasSourceTerms() const;
        
        /*
         * Set the weak pointer to the flow model from the parent FlowModel class.
         */
        void setFlowModel(const HAMERS_WEAK_PTR<FlowModel>& flow_model)
        {
            d_flow_model = flow_model;
        }
        
        /*
         * Register the required variables for the computation of source terms in the registered patch.
         */
        virtual void
        registerDerivedVariablesForSourceTerms(
            const hier::IntVector& num_subghosts);
        
        /*
         * Register the required variables for the computation of local stable time increment for
         * source terms in the registered patch.
         */
        virtual void
        registerDerivedVariablesForSourceTermsStableDt(
            const hier::IntVector& num_subghosts);
        
        /*
         * Allocate memory for cell data of different registered derived variables related to this
         * class in the registered patch.
         */
        virtual void allocateMemoryForDerivedCellData();
        
        /*
         * Clear cell data of different derived variables related to this class in the registered patch.
         */
        virtual void clearCellData();
        
        /*
         * Compute cell data of different registered derived variables related to this class.
         */
        virtual void computeDerivedCellData();
        
        /*
         * Compute all source terms on a patch.
         */
        virtual void
        computeSourceTermsOnPatch(
            const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_source,
            const double time,
            const double dt,
            const int RK_step_number);
        
        /*
         * Get local stable time increment for source terms.
         */
        virtual double
        getStableDtOnPatch();
        
        /*
         * Put the characteristics of this class into the restart database.
         */
        virtual void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
    protected:
        /*
         * Compute special source terms.
         */
        void
        computeSpecialSourceTermsOnPatch(
            const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_source,
            const double time,
            const double dt,
            const int RK_step_number);
        
        /*
         * Put the characteristics of base class into the restart database.
         */
        void
        putToRestartBase(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Put the characteristics of base class into the restart source database.
         */
        void
        putToRestartSourceBase(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_source_terms_db) const;
        
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Name of the project.
         */
        std::string d_project_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * Whether there are source terms.
         */
        bool d_has_source_terms;
        
        /*
         * Whether there are special source terms.
         */
        bool d_has_special_source_terms;
        
        /*
         * Special source terms object.
         */
        HAMERS_SHARED_PTR<FlowModelSpecialSourceTerms> d_special_source_terms;
        
        /* 
         * Whether all derived cell data related to this class is computed in full domain or sub-domain.
         */
        bool d_derived_cell_data_computed;
        
        /*
         * Number of sub-ghost cells of source terms.
         */
        hier::IntVector d_num_subghosts_source_terms;
        
        /*
         * Box with sub-ghost cells of source terms.
         */
        hier::Box d_subghost_box_source_terms;
        
        /*
         * Dimensions of box with sub-ghost cells of source terms.
         */
        hier::IntVector d_subghostcell_dims_source_terms;
        
        /*
         * HAMERS_WEAK_PTR to FlowModel.
         */
        HAMERS_WEAK_PTR<FlowModel> d_flow_model;
        
};

#endif /* FLOW_MODEL_SOURCE_UTILITIES_HPP */
