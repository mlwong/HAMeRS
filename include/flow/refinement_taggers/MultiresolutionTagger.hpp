#ifndef MULTIRESOLUTION_TAGGER_HPP
#define MULTIRESOLUTION_TAGGER_HPP

#include "HAMeRS_config.hpp"

#include "algs/integrator/RungeKuttaLevelIntegrator.hpp"
#include "flow/flow_models/FlowModels.hpp"
#include "util/wavelet_transform/WaveletTransformHarten.hpp"

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
         * Compute values of multiresolution sensors on a patch.
         */
        void
        computeMultiresolutionSensorValuesOnPatch(
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
         * Tag cells on a patch for refinement using multiresolution sensors.
         */
        void
        tagCellsOnPatch(
            hier::Patch& patch,
            const boost::shared_ptr<pdat::CellData<int> >& tags,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
    private:
        /*
         * Tag cells on a patch using wavelet sensor with the combination of three possible criteria:
         * 1. When ratio between wavelet coefficient and global maximum at any level is greater than the tolerance.
         * 2. When ratio between wavelet coefficient and local mean at any level is greater than the tolerance.
         * 3. When the Lipschitz's exponent is smaller than the tolerance.
         */
        void
        tagCellsOnPatchWithWaveletSensor(
            hier::Patch& patch,
            const boost::shared_ptr<pdat::CellData<int> >& tags,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& wavelet_coeffs,
            const std::vector<double>& wavelet_coeffs_maxs,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& variable_local_means,
            const boost::shared_ptr<pdat::CellData<double> >& Lipschitz_exponent,
            const std::string& sensor_key,
            const bool uses_global_tol,
            const bool uses_local_tol,
            const bool uses_alpha_tol,
            const double global_tol,
            const double local_tol,
            const double alpha_tol);
        
        /*
         * Compute the Lipschitz's exponent on a patch. There are two steps:
         * 1. Find the maximum wavelet coefficients in the domain of dependence.
         * 2. Compute Lipschitz's exponent.
         */
        void
        computeLipschitzExponentOnPatch(
            hier::Patch& patch,
            const boost::shared_ptr<pdat::CellData<double> >& Lipschitz_exponent,
            const std::vector<boost::shared_ptr<pdat::CellData<double> > >& wavelet_coeffs,
            const std::string& sensor_key);
        
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
         * Number of ghost cells needed by the the multiresolution tagger.
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
         * Number of levels and number of vanishing moments for WaveletTransformHarten.
         */
        int d_Harten_wavelet_num_level;
        int d_Harten_wavelet_num_vanishing_moments;
        
        /*
         * Variables, tolerances and settings for the multiresolution sensors.
         */
        std::vector<std::string> d_Harten_wavelet_variables;
        std::vector<double> d_Harten_wavelet_global_tol;
        std::vector<double> d_Harten_wavelet_local_tol;
        std::vector<double> d_Harten_wavelet_alpha_tol;
        std::vector<bool> d_Harten_wavelet_uses_global_tol;
        std::vector<bool> d_Harten_wavelet_uses_local_tol;
        std::vector<bool> d_Harten_wavelet_uses_alpha_tol;
        
        /*
         * boost::shared_ptr to wavelet coefficients at different levels.
         */
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_wavelet_coeffs_density;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_wavelet_coeffs_total_energy;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_wavelet_coeffs_pressure;
        
        /*
         * Statistics of sensor values.
         */
        std::vector<double> d_Harten_wavelet_coeffs_maxs_density;
        std::vector<double> d_Harten_wavelet_coeffs_maxs_total_energy;
        std::vector<double> d_Harten_wavelet_coeffs_maxs_pressure;
        
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_local_means_density;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_local_means_total_energy;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_local_means_pressure;
        
        boost::shared_ptr<pdat::CellVariable<double> > d_Harten_Lipschitz_exponent_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_Harten_Lipschitz_exponent_total_energy;
        boost::shared_ptr<pdat::CellVariable<double> > d_Harten_Lipschitz_exponent_pressure;
        
};

#endif /* MULTIRESOLUTION_TAGGER_HPP */
