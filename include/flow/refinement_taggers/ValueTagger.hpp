#ifndef VALUE_TAGGER_HPP
#define VALUE_TAGGER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "algs/integrator/RungeKuttaLevelIntegrator.hpp"
#include "extn/visit_data_writer/ExtendedVisItDataWriter.hpp"
#include "flow/flow_models/FlowModels.hpp"
#include "util/derivatives/DerivativeFirstOrder.hpp"

// #include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

class ValueTagger
{
    public:
        ValueTagger(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& value_tagger_db);
        
        /*
         * Get the number of ghost cells needed by the value tagger.
         */
        hier::IntVector
        getValueTaggerNumberOfGhostCells() const
        {
            return d_num_value_ghosts;
        }
        
        /*
         * Register the variables used in value tagger class.
         */
        void
        registerValueTaggerVariables(
            RungeKuttaLevelIntegrator* integrator);
        
        /*
         * Register the plotting quantities.
         */
        void
        registerPlotQuantities(
            const HAMERS_SHARED_PTR<ExtendedVisItDataWriter>& visit_writer,
            const HAMERS_SHARED_PTR<hier::VariableContext>& plot_context);
        
        /*
         * Print all characteristics of the value tagger class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the value tagger class into the restart
         * database.
         */
        void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Compute values on a patch for value tagger.
         */
        void
        computeValueTaggerValuesOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Get the statistics of values that are required by the value tagger.
         */
        void
        getValueStatistics(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Tag cells on a patch for refinement using value tagger.
         */
        void
        tagCellsOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
    private:
        /*
         * Tag cells on a patch for refinement using data values.
         */
        void
        tagCellsOnPatchWithValue(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags,
            const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_value_tagger,
            const Real value_max,
            const bool uses_global_tol_up,
            const bool uses_global_tol_lo,
            const bool uses_local_tol_up,
            const bool uses_local_tol_lo,
            const Real global_tol_up,
            const Real global_tol_lo,
            const Real local_tol_up,
            const Real local_tol_lo);
        
        /*
         * Transfer data input on a patch to data in class variable.
         */
        void transferDataOnPatchToClassVariable(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& data_input,
            const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_value_tagger,
            const int depth);
        
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * HAMERS_SHARED_PTR to the grid geometry.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of ghost cells needed by the the value detector.
         */
        hier::IntVector d_num_value_ghosts;
        
        /*
         * Flow model.
         */
        const HAMERS_SHARED_PTR<FlowModel> d_flow_model;
        
        /*
         * Number of ghost cells to use in taking derivatives.
         */
        int d_num_ghosts_derivative;
        
        /*
         * Variables, tolerances and settings for the value sensor.
         */
        std::vector<std::string> d_variables;
        std::vector<Real> d_global_tol_up;
        std::vector<Real> d_global_tol_lo;
        std::vector<Real> d_local_tol_up;
        std::vector<Real> d_local_tol_lo;
        std::vector<bool> d_uses_global_tol_up;
        std::vector<bool> d_uses_global_tol_lo;
        std::vector<bool> d_uses_local_tol_up;
        std::vector<bool> d_uses_local_tol_lo;
        
        /*
         * HAMERS_SHARED_PTR to data values.
         */
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_value_tagger_variable_density;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_value_tagger_variable_total_energy;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_value_tagger_variable_pressure;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_value_tagger_variable_dilatation;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_value_tagger_variable_enstrophy;
        std::vector<HAMERS_SHARED_PTR<pdat::CellVariable<Real> > > d_value_tagger_variable_mass_fractions;
        
        /*
         * Statistics of data values.
         */
        Real d_value_tagger_max_density;
        Real d_value_tagger_max_total_energy;
        Real d_value_tagger_max_pressure;
        Real d_value_tagger_max_dilatation;
        Real d_value_tagger_max_enstrophy;
        std::vector<Real> d_value_tagger_max_mass_fractions;
        
};

#endif /* VALUE_TAGGER_HPP */
