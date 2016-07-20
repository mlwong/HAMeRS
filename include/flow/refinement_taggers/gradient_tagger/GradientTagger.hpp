#ifndef GRADIENT_TAGGER_HPP
#define GRADIENT_TAGGER_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/VisItDataWriter.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Dimension.h"

#include "algs/integrator/RungeKuttaLevelIntegrator.hpp"
#include "flow/flow_models/FlowModels.hpp"
#include "util/gradient_sensors/GradientSensorJameson.hpp"

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
         * Tag cells for refinement using gradient sensors.
         */
        void
        tagCells(
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<int> > tags,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
    private:
        /*
         * Tag cells using gradient sensor.
         */
        void
        tagCellsWithGradientSensor(
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<int> > tags,
            boost::shared_ptr<pdat::CellData<double> > gradient,
            double& tol,
            std::string& sensor_key);
        
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
         * Number of ghost cells needed by the the gradient detector.
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
         * boost::shared_ptr to GradientSensorJameson.
         */
        boost::shared_ptr<GradientSensorJameson> d_gradient_sensor_Jameson;
        
        /*
         * Variables and tolerances for the gradient sensor.
         */
        std::vector<std::string> d_Jameson_gradient_variables;
        std::vector<double> d_Jameson_gradient_tol;
        
        /*
         * boost::shared_ptr to Jameson's gradient.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_Jameson_density_gradient;
        boost::shared_ptr<pdat::CellVariable<double> > d_Jameson_total_energy_gradient;
        boost::shared_ptr<pdat::CellVariable<double> > d_Jameson_pressure_gradient;
        boost::shared_ptr<pdat::CellVariable<double> > d_Jameson_dilatation_gradient;
        boost::shared_ptr<pdat::CellVariable<double> > d_Jameson_enstrophy_gradient;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Jameson_mass_fraction_gradient;
        
};

#endif /* GRADIENT_TAGGER_HPP */
