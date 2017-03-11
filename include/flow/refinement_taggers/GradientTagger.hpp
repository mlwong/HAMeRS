#ifndef GRADIENT_TAGGER_HPP
#define GRADIENT_TAGGER_HPP

#include "HAMeRS_config.hpp"

#include "algs/integrator/RungeKuttaLevelIntegrator.hpp"
#include "flow/flow_models/FlowModels.hpp"
#include "util/differences/DifferenceFirstDerivative.hpp"
#include "util/differences/DifferenceSecondDerivative.hpp"
#include "util/gradient_sensors/GradientSensorJameson.hpp"

#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Dimension.h"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

class GradientTagger
{
    public:
        GradientTagger(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& gradient_tagger_db);
        
        /*
         * Get the number of ghost cells needed by the gradient tagger.
         */
        hier::IntVector
        getGradientTaggerNumberOfGhostCells() const
        {
            return d_num_gradient_ghosts;
        }
        
        /*
         * Register the variables used in gradient tagger class.
         */
        void
        registerGradientTaggerVariables(
            RungeKuttaLevelIntegrator* integrator);
        
        /*
         * Register the plotting quantities.
         */
        void
        registerPlotQuantities(
            const boost::shared_ptr<appu::VisItDataWriter>& visit_writer,
            const boost::shared_ptr<hier::VariableContext>& plot_context);
        
        /*
         * Print all characteristics of the gradient tagger class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the gradient tagger class into the restart
         * database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Compute values of gradient sensors on patch.
         */
        void
        computeGradientSensorValuesOnPatch(
            hier::Patch& patch,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Get the statistics of the sensor values that are required by the
         * gradient sensors at a given patch level.
         */
        void
        getSensorValueStatistics(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Tag cells on patch for refinement using gradient sensors.
         */
        void
        tagCellsOnPatch(
            hier::Patch& patch,
            const boost::shared_ptr<pdat::CellData<int> >& tags,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
    private:
        /*
         * Tag cells on patch using value of gradient sensor.
         */
        void
        tagCellsOnPatchWithGradientSensor(
            hier::Patch& patch,
            const boost::shared_ptr<pdat::CellData<int> >& tags,
            const boost::shared_ptr<pdat::CellData<double> >& gradient,
            const std::string& sensor_key,
            const double tol);
        
        /*
         * Tag cells on patch using difference sensor.
         */
        void
        tagCellsOnPatchWithDifferenceSensor(
            hier::Patch& patch,
            const boost::shared_ptr<pdat::CellData<int> >& tags,
            const boost::shared_ptr<pdat::CellData<double> >& difference,
            const double difference_max,
            const boost::shared_ptr<pdat::CellData<double> >& variable_local_mean,
            const bool uses_global_tol,
            const bool uses_local_tol,
            const double global_tol,
            const double local_tol);
        
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
         * Number of ghost cells needed by the the gradient tagger.
         */
        hier::IntVector d_num_gradient_ghosts;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Flow model.
         */
        const boost::shared_ptr<FlowModel> d_flow_model;
        
        /*
         * Chosen gradient sensors.
         */
        std::vector<std::string> d_gradient_sensors;
        
        /*
         * boost::shared_ptr to difference operators.
         */
        boost::shared_ptr<DifferenceFirstDerivative> d_difference_first_derivative;
        boost::shared_ptr<DifferenceSecondDerivative> d_difference_second_derivative;
        
        /*
         * boost::shared_ptr to GradientSensorJameson.
         */
        boost::shared_ptr<GradientSensorJameson> d_gradient_sensor_Jameson;
        
        /*
         * Variables, tolerances and settings for the gradient sensors.
         */
        std::vector<std::string> d_difference_first_derivative_variables;
        std::vector<double> d_difference_first_derivative_global_tol;
        std::vector<double> d_difference_first_derivative_local_tol;
        std::vector<bool> d_difference_first_derivative_uses_global_tol;
        std::vector<bool> d_difference_first_derivative_uses_local_tol;
        
        std::vector<std::string> d_difference_second_derivative_variables;
        std::vector<double> d_difference_second_derivative_global_tol;
        std::vector<double> d_difference_second_derivative_local_tol;
        std::vector<bool> d_difference_second_derivative_uses_global_tol;
        std::vector<bool> d_difference_second_derivative_uses_local_tol;
        
        /*
         * Variables and tolerances for the Jameson gradient sensors.
         */
        std::vector<std::string> d_Jameson_gradient_variables;
        std::vector<double> d_Jameson_gradient_tol;
        
        /*
         * boost::shared_ptr to differences.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_difference_first_derivative_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_difference_first_derivative_total_energy;
        boost::shared_ptr<pdat::CellVariable<double> > d_difference_first_derivative_pressure;
        
        boost::shared_ptr<pdat::CellVariable<double> > d_difference_second_derivative_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_difference_second_derivative_total_energy;
        boost::shared_ptr<pdat::CellVariable<double> > d_difference_second_derivative_pressure;
        
        /*
         * boost::shared_ptr to values of Jameson gradient sensor.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_Jameson_gradient_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_Jameson_gradient_total_energy;
        boost::shared_ptr<pdat::CellVariable<double> > d_Jameson_gradient_pressure;
        
        /*
         * Statistics of sensor values.
         */
        double d_difference_first_derivative_max_density;
        double d_difference_first_derivative_max_total_energy;
        double d_difference_first_derivative_max_pressure;
        
        double d_difference_second_derivative_max_density;
        double d_difference_second_derivative_max_total_energy;
        double d_difference_second_derivative_max_pressure;
        
        boost::shared_ptr<pdat::CellVariable<double> > d_difference_first_derivative_local_mean_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_difference_first_derivative_local_mean_total_energy;
        boost::shared_ptr<pdat::CellVariable<double> > d_difference_first_derivative_local_mean_pressure;
        
        boost::shared_ptr<pdat::CellVariable<double> > d_difference_second_derivative_local_mean_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_difference_second_derivative_local_mean_total_energy;
        boost::shared_ptr<pdat::CellVariable<double> > d_difference_second_derivative_local_mean_pressure;
        
};

#endif /* GRADIENT_TAGGER_HPP */
