#ifndef VALUE_TAGGER_HPP
#define VALUE_TAGGER_HPP

#include "HAMeRS_config.hpp"

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

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

class ValueTagger
{
    public:
        ValueTagger(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& value_tagger_db);
        
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
            const boost::shared_ptr<ExtendedVisItDataWriter>& visit_writer,
            const boost::shared_ptr<hier::VariableContext>& plot_context);
        
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
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Compute values on a patch for value tagger.
         */
        void
        computeValueTaggerValuesOnPatch(
            hier::Patch& patch,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Get the statistics of values that are required by the value tagger.
         */
        void
        getValueStatistics(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Tag cells on a patch for refinement using value tagger.
         */
        void
        tagCellsOnPatch(
            hier::Patch& patch,
            const boost::shared_ptr<pdat::CellData<int> >& tags,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
    private:
        /*
         * Tag cells on a patch for refinement using data values.
         */
        void
        tagCellsOnPatchWithValue(
            hier::Patch& patch,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const boost::shared_ptr<pdat::CellData<int> >& tags,
            const boost::shared_ptr<pdat::CellVariable<double> >& variable_value_tagger,
            const double value_max,
            const bool uses_global_tol_up,
            const bool uses_global_tol_lo,
            const bool uses_local_tol_up,
            const bool uses_local_tol_lo,
            const double global_tol_up,
            const double global_tol_lo,
            const double local_tol_up,
            const double local_tol_lo);
        
        /*
         * Transfer data input on a patch to data in class variable.
         */
        void transferDataOnPatchToClassVariable(
            hier::Patch& patch,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const boost::shared_ptr<pdat::CellData<double> >& data_input,
            const boost::shared_ptr<pdat::CellVariable<double> >& variable_value_tagger,
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
         * boost::shared_ptr to the grid geometry.
         */
        const boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Number of ghost cells needed by the the value detector.
         */
        hier::IntVector d_num_value_ghosts;
        
        /*
         * Flow model.
         */
        const boost::shared_ptr<FlowModel> d_flow_model;
        
        /*
         * Number of ghost cells to use in taking derivatives.
         */
        int d_num_ghosts_derivative;
        
        /*
         * Variables, tolerances and settings for the value sensor.
         */
        std::vector<std::string> d_variables;
        std::vector<double> d_global_tol_up;
        std::vector<double> d_global_tol_lo;
        std::vector<double> d_local_tol_up;
        std::vector<double> d_local_tol_lo;
        std::vector<bool> d_uses_global_tol_up;
        std::vector<bool> d_uses_global_tol_lo;
        std::vector<bool> d_uses_local_tol_up;
        std::vector<bool> d_uses_local_tol_lo;
        
        /*
         * boost::shared_ptr to data values.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_value_tagger_variable_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_value_tagger_variable_total_energy;
        boost::shared_ptr<pdat::CellVariable<double> > d_value_tagger_variable_pressure;
        boost::shared_ptr<pdat::CellVariable<double> > d_value_tagger_variable_dilatation;
        boost::shared_ptr<pdat::CellVariable<double> > d_value_tagger_variable_enstrophy;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_value_tagger_variable_mass_fractions;
        
        /*
         * Statistics of data values.
         */
        double d_value_tagger_max_density;
        double d_value_tagger_max_total_energy;
        double d_value_tagger_max_pressure;
        double d_value_tagger_max_dilatation;
        double d_value_tagger_max_enstrophy;
        std::vector<double> d_value_tagger_max_mass_fractions;
        
};

#endif /* VALUE_TAGGER_HPP */
