#ifndef GRADIENT_TAGGER_HPP
#define GRADIENT_TAGGER_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "algs/integrator/RungeKuttaLevelIntegrator.hpp"
#include "extn/visit_data_writer/ExtendedVisItDataWriter.hpp"
#include "flow/flow_models/FlowModels.hpp"
#include "util/differences/DifferenceFirstOrder.hpp"
#include "util/differences/DifferenceSecondOrder.hpp"
#include "util/gradient_sensors/GradientSensorJameson.hpp"
#include "util/gradient_sensors/GradientSensorDucros.hpp"

// #include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <vector>

class GradientTagger
{
    public:
        GradientTagger(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& gradient_tagger_db);
        
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
            const HAMERS_SHARED_PTR<ExtendedVisItDataWriter>& visit_writer,
            const HAMERS_SHARED_PTR<hier::VariableContext>& plot_context);
        
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
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Compute values of gradient sensors on a patch.
         */
        void
        computeGradientSensorValuesOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Get the statistics of the sensor values that are required by the
         * gradient sensors at a given patch level.
         */
        void
        getSensorValueStatistics(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
        /*
         * Tag cells on a patch for refinement using gradient sensors.
         */
        void
        tagCellsOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context);
        
    private:
        /*
         * Tag cells on a patch using value of gradient sensor.
         */
        void
        tagCellsOnPatchWithGradientSensor(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& gradient,
            const std::string& sensor_key,
            const Real tol);
        
        /*
         * Tag cells on a patch using difference sensor.
         */
        void
        tagCellsOnPatchWithDifferenceSensor(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::CellData<int> >& tags,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& difference,
            const Real difference_max,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& variable_local_mean,
            const bool uses_global_tol,
            const bool uses_local_tol,
            const Real global_tol,
            const Real local_tol);
        
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
         * Number of ghost cells needed by the the gradient tagger.
         */
        hier::IntVector d_num_gradient_ghosts;
        
        /*
         * Flow model.
         */
        const HAMERS_SHARED_PTR<FlowModel> d_flow_model;
        
        /*
         * Chosen gradient sensors.
         */
        std::vector<std::string> d_gradient_sensors;
        
        /*
         * Whether only apply the gradient sensors in a box.
         */
        bool d_is_tagging_in_box_only;
        std::vector<Real> d_tagging_box_xlo; // Lower spatial coordinates.
        std::vector<Real> d_tagging_box_xhi; // Upper spatial coordinates.
        
        /*
         * HAMERS_SHARED_PTR to difference operators.
         */
        HAMERS_SHARED_PTR<DifferenceFirstOrder> d_difference_first_order;
        HAMERS_SHARED_PTR<DifferenceSecondOrder> d_difference_second_order;
        
        /*
         * HAMERS_SHARED_PTR to GradientSensorJameson.
         */
        HAMERS_SHARED_PTR<GradientSensorJameson> d_gradient_sensor_Jameson;
        
        /*
         * HAMERS_SHARED_PTR to GradientSensorDucros.
         */
        HAMERS_SHARED_PTR<GradientSensorDucros> d_gradient_sensor_Ducros;
        
        /*
         * Variables, tolerances and settings for the gradient sensors.
         */
        std::vector<std::string> d_difference_first_order_variables;
        std::vector<Real> d_difference_first_order_global_tol;
        std::vector<Real> d_difference_first_order_local_tol;
        std::vector<bool> d_difference_first_order_uses_global_tol;
        std::vector<bool> d_difference_first_order_uses_local_tol;
        
        std::vector<std::string> d_difference_second_order_variables;
        std::vector<Real> d_difference_second_order_global_tol;
        std::vector<Real> d_difference_second_order_local_tol;
        std::vector<bool> d_difference_second_order_uses_global_tol;
        std::vector<bool> d_difference_second_order_uses_local_tol;
        
        /*
         * Variables and tolerances for the Jameson gradient sensors.
         */
        std::vector<std::string> d_Jameson_gradient_variables;
        std::vector<Real> d_Jameson_gradient_tol;
        
        /*
         * Settings and tolerance for the Ducros gradient sensor.
         */
        bool d_Ducros_gradient_use_strain_rate_instead_of_dilatation;
        Real d_Ducros_gradient_tol;
        
        /*
         * HAMERS_SHARED_PTR to differences.
         */
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_difference_first_order_density;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_difference_first_order_total_energy;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_difference_first_order_pressure;
        
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_difference_second_order_density;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_difference_second_order_total_energy;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_difference_second_order_pressure;
        
        /*
         * HAMERS_SHARED_PTR to values of Jameson gradient sensor.
         */
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_Jameson_gradient_density;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_Jameson_gradient_total_energy;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_Jameson_gradient_pressure;
        
        /*
         * HAMERS_SHARED_PTR to of value of Ducros gradient sensor.
         */
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_Ducros_gradient;
         
        /*
         * Statistics of sensor values.
         */
        Real d_difference_first_order_max_density;
        Real d_difference_first_order_max_total_energy;
        Real d_difference_first_order_max_pressure;
        
        Real d_difference_second_order_max_density;
        Real d_difference_second_order_max_total_energy;
        Real d_difference_second_order_max_pressure;
        
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_difference_first_order_local_mean_density;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_difference_first_order_local_mean_total_energy;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_difference_first_order_local_mean_pressure;
        
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_difference_second_order_local_mean_density;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_difference_second_order_local_mean_total_energy;
        HAMERS_SHARED_PTR<pdat::CellVariable<Real> > d_difference_second_order_local_mean_pressure;
        
};

#endif /* GRADIENT_TAGGER_HPP */
