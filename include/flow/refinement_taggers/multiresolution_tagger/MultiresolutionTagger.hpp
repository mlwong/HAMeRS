#ifndef MULTIRESOLUTION_TAGGER_HPP
#define MULTIRESOLUTION_TAGGER_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/VisItDataWriter.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Dimension.h"

#include "algs/integrator/RungeKuttaLevelIntegrator.hpp"
#include "flow/flow_models/FlowModels.hpp"
#include "util/wavelet_transform/WaveletTransformHarten.hpp"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

class MultiresolutionTagger
{
    public:
        MultiresolutionTagger(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_species,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& multiresolution_tagger_db);
        
        /*
         * Get the number of ghost cells needed by the multiresolution tagger.
         */
        hier::IntVector
        getMultiresolutionTaggerNumberOfGhostCells() const
        {
            return d_num_multiresolution_ghosts;
        }
        
        /*
         * Register the variables used in multiresolution tagger class.
         */
        void
        registerMultiresolutionTaggerVariables(
            RungeKuttaLevelIntegrator* integrator);
        
        /*
         * Register the plotting quantities.
         */
        void
        registerPlotQuantities(
            const boost::shared_ptr<appu::VisItDataWriter>& visit_writer,
            const boost::shared_ptr<hier::VariableContext>& plot_context);
        
        /*
         * Print all characteristics of the multiresolution tagger class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the multiresolution tagger into the restart
         * database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Compute values of multiresolution sensors.
         */
        void
        computeMultiresolutionSensorValues(
            hier::Patch& patch,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Get the statistics of the sensor values that are required by the
         * multiresolution sensors at a given patch level.
         */
        void
        getSensorValueStatistics(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Tag cells for refinement using multiresolution sensors.
         */
        void
        tagCells(
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<int> > tags,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
    private:
        /*
         * Tag cells using wavelet sensor with the combination of three possible criteria:
         * 1. When ratio between wavelet coefficient and global maximum at any level is greater than the tolerance.
         * 2. When ratio between wavelet coefficient and local mean at any level is greater than the tolerance.
         * 3. When the Lipschitz's exponent is smaller than the tolerance.
         */
        void
        tagCellsWithWaveletSensor(
            hier::Patch& patch,
            boost::shared_ptr<pdat::CellData<int> > tags,
            std::vector<boost::shared_ptr<pdat::CellData<double> > >& wavelet_coeffs,
            std::vector<double>& wavelet_coeffs_maxs,
            std::vector<boost::shared_ptr<pdat::CellData<double> > >& variable_local_means,
            boost::shared_ptr<pdat::CellData<double> > Lipschitz_exponent,
            double& global_tol,
            double& local_tol,
            double& alpha_tol,
            std::string& sensor_key);
        
        /*
         * Compute the Lipschitz's exponent. There are two steps:
         * 1. Find the maximum wavelet coefficients in the domain of dependence.
         * 2. Compute Lipschitz's exponent.
         */
        void
        computeLipschitzExponent(
            hier::Patch& patch,
            std::vector<boost::shared_ptr<pdat::CellData<double> > >& wavelet_coeffs,
            boost::shared_ptr<pdat::CellData<double> > Lipschitz_exponent,
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
         * Number of ghost cells needed by the the multiresolution detector.
         */
        hier::IntVector d_num_multiresolution_ghosts;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * Flow model.
         */
        const boost::shared_ptr<FlowModel> d_flow_model;
        
        /*
         * Chosen multiresolution sensors.
         */
        std::vector<std::string> d_multiresolution_sensors;
        
        /*
         * boost::shared_ptr to WaveletTransformHarten.
         */
        boost::shared_ptr<WaveletTransformHarten> d_wavelet_transfrom_Harten;
        
        /*
         * Number of levels and number of vanishing moments of WaveletTransformHarten.
         */
        int d_Harten_wavelet_num_level;
        int d_Harten_wavelet_num_vanishing_moments;
        
        /*
         * Variables, tolerances and settings for the multiresolution sensor.
         */
        std::vector<std::string> d_Harten_wavelet_variables;
        std::vector<double> d_Harten_wavelet_global_tol;
        std::vector<double> d_Harten_wavelet_local_tol;
        std::vector<double> d_Harten_wavelet_alpha_tol;
        bool d_Harten_wavelet_uses_global_tol;
        bool d_Harten_wavelet_uses_local_tol;
        bool d_Harten_wavelet_uses_alpha_tol;
        
        /*
         * boost::shared_ptr to wavelet coefficients at different levels.
         */
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_density_wavelet_coeffs;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_total_energy_wavelet_coeffs;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_pressure_wavelet_coeffs;
        
        /*
         * Statistics of sensor values.
         */
        std::vector<double> d_Harten_density_wavelet_coeffs_maxs;
        std::vector<double> d_Harten_total_energy_wavelet_coeffs_maxs;
        std::vector<double> d_Harten_pressure_wavelet_coeffs_maxs;
        
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_density_local_means;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_total_energy_local_means;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_pressure_local_means;
        
        boost::shared_ptr<pdat::CellVariable<double> > d_Harten_density_Lipschitz_exponent;
        boost::shared_ptr<pdat::CellVariable<double> > d_Harten_total_energy_Lipschitz_exponent;
        boost::shared_ptr<pdat::CellVariable<double> > d_Harten_pressure_Lipschitz_exponent;
        
};

#endif /* MULTIRESOLUTION_TAGGER_HPP */
