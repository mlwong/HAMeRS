#ifndef FLOW_MODEL_SOURCE_UTILITIES_HPP
#define FLOW_MODEL_SOURCE_UTILITIES_HPP

#include "HAMeRS_config.hpp"

#include "flow/flow_models/FlowModel.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellData.h"

#include "boost/weak_ptr.hpp"
#include <string>

class FlowModel;

class FlowModelSourceUtilities
{
    public:
        FlowModelSourceUtilities(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const int& num_eqn,
            const boost::shared_ptr<tbox::Database>& flow_model_db);
        
        virtual ~FlowModelSourceUtilities() {}
        
        /*
         * Check whether there are any source terms.
         */
        bool hasSourceTerms() const;
        
        /*
         * Set the weak pointer to the flow model from the parent FlowModel class.
         */
        void setFlowModel(const boost::weak_ptr<FlowModel>& flow_model)
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
            const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
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
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
protected:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * boost::shared_ptr to the grid geometry.
         */
        const boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
        
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
         * boost::weak_ptr to FlowModel.
         */
        boost::weak_ptr<FlowModel> d_flow_model;
        
};

#endif /* FLOW_MODEL_SOURCE_UTILITIES_HPP */
