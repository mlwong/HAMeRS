#ifndef MULTIRESOLUTION_TAGGER_HPP
#define MULTIRESOLUTION_TAGGER_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/VisItDataWriter.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Dimension.h"

#include "equation_of_state/EquationOfStateIdealGas.hpp"
#include "flow_model/FlowModels.hpp"
#include "integrator/RungeKuttaLevelIntegrator.hpp"
#include "utils/wavelet_transform/WaveletTransformHarten.hpp"

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
            const FLOW_MODEL& flow_model,
            const int& num_species,
            const boost::shared_ptr<EquationOfState>& equation_of_state,
            const boost::shared_ptr<tbox::Database>& multiresolution_tagger_db);
        
        /*
         * Get the number of ghost cells needed by the multiresolution tagger.
         */
        hier::IntVector
        getMultiresolutionTaggerNumberOfGhostCells() const
        {
            hier::IntVector d_num_multiresolution_ghosts = hier::IntVector::getZero(d_dim);
            
            if (d_wavelet_transfrom_Harten != nullptr)
            {
                d_num_multiresolution_ghosts = hier::IntVector::max(
                    d_num_multiresolution_ghosts,
                    d_wavelet_transfrom_Harten->getWaveletTransformNumberOfGhostCells());
            }
            
            return d_num_multiresolution_ghosts;
        }
        
        /*
         * Set the number of ghost cells needed.
         */
        void
        setNumberOfGhostCells(const hier::IntVector& num_ghosts)
        {
            d_num_ghosts = num_ghosts;
            
            if (d_wavelet_transfrom_Harten != nullptr)
            {
                d_wavelet_transfrom_Harten->setNumberOfGhostCells(num_ghosts);
            }
            
            d_num_ghosts_set = true;
        }
        
        /*
         * Set the cell variables if single-species flow model is chosen.
         */
        void
        setVariablesForSingleSpecies(
            const boost::shared_ptr<pdat::CellVariable<double> >& density,
            const boost::shared_ptr<pdat::CellVariable<double> >& momentum,
            const boost::shared_ptr<pdat::CellVariable<double> >& total_energy)
        {
            if (d_flow_model != SINGLE_SPECIES)
            {            
                TBOX_ERROR(d_object_name
                           << ": "
                           << "setVariablesForSingleSpecies() shouldn't be used."
                           << std::endl);
            }
            
            d_density = density;
            d_momentum = momentum;
            d_total_energy = total_energy;
            
            d_variables_set = true;
        }
        
        /*
         * Set the cell variables if four-equation multi-species flow model
         * by Shyue is chosen.
         */
        void
        setVariablesForFourEqnShyue(
            const boost::shared_ptr<pdat::CellVariable<double> >& density,
            const boost::shared_ptr<pdat::CellVariable<double> >& momentum,
            const boost::shared_ptr<pdat::CellVariable<double> >& total_energy,
            const boost::shared_ptr<pdat::CellVariable<double> >& mass_fraction)
        {
            if (d_flow_model != FOUR_EQN_SHYUE)
            {            
                TBOX_ERROR(d_object_name
                           << ": "
                           << "setVariablesForFourEqnShyue() shouldn't be used."
                           << std::endl);
            }
            
            d_density = density;
            d_momentum = momentum;
            d_total_energy = total_energy;
            d_mass_fraction = mass_fraction;
            
            d_variables_set = true;
        }
        
        /*
         * Set the cell variables if five-equation multi-species flow model
         * by Allaire is chosen.
         */
        void
        setVariablesForFiveEqnAllaire(
            const boost::shared_ptr<pdat::CellVariable<double> >& partial_density,
            const boost::shared_ptr<pdat::CellVariable<double> >& momentum,
            const boost::shared_ptr<pdat::CellVariable<double> >& total_energy,
            const boost::shared_ptr<pdat::CellVariable<double> >& volume_fraction)
        {
            if (d_flow_model != FIVE_EQN_ALLAIRE)
            {            
                TBOX_ERROR(d_object_name
                           << ": "
                           << "setVariablesForFiveEqnAllaire() shouldn't be used."
                           << std::endl);
            }
            
            d_partial_density = partial_density;
            d_momentum = momentum;
            d_total_energy = total_energy;
            d_volume_fraction = volume_fraction;
            
            d_variables_set = true;
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
            RungeKuttaLevelIntegrator* integrator,
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
         * Number of ghost cells for time-independent variables.
         */
        hier::IntVector d_num_ghosts;
        
        /*
         * Flow model.
         */
        const FLOW_MODEL d_flow_model;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * boost::shared_ptr to EquationOfState.
         */
        const boost::shared_ptr<EquationOfState> d_equation_of_state;
        
        /*
         * boost::shared_ptr to solution state.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_partial_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_momentum;
        boost::shared_ptr<pdat::CellVariable<double> > d_total_energy;
        boost::shared_ptr<pdat::CellVariable<double> > d_mass_fraction;
        boost::shared_ptr<pdat::CellVariable<double> > d_volume_fraction;
        
        /*
         * Boolean to determine whether proper variables are initialized.
         */
        bool d_variables_set;
        
        /*
         * Boolean to determine whether the number of ghost cells is initialized.
         */
        bool d_num_ghosts_set;
        
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
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_enstrophy_wavelet_coeffs;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_mass_fraction_wavelet_coeffs;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_volume_fraction_wavelet_coeffs;
        
        /*
         * Statistics of sensor values.
         */
        std::vector<double> d_Harten_density_wavelet_coeffs_maxs;
        std::vector<double> d_Harten_total_energy_wavelet_coeffs_maxs;
        std::vector<double> d_Harten_pressure_wavelet_coeffs_maxs;
        std::vector<double> d_Harten_enstrophy_wavelet_coeffs_maxs;
        std::vector<double> d_Harten_mass_fraction_wavelet_coeffs_maxs;
        std::vector<double> d_Harten_volume_fraction_wavelet_coeffs_maxs;
        
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_density_local_means;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_total_energy_local_means;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_pressure_local_means;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_enstrophy_local_means;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_mass_fraction_local_means;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_volume_fraction_local_means;
        
        boost::shared_ptr<pdat::CellVariable<double> > d_Harten_density_Lipschitz_exponent;
        boost::shared_ptr<pdat::CellVariable<double> > d_Harten_total_energy_Lipschitz_exponent;
        boost::shared_ptr<pdat::CellVariable<double> > d_Harten_pressure_Lipschitz_exponent;
        boost::shared_ptr<pdat::CellVariable<double> > d_Harten_enstrophy_Lipschitz_exponent;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_mass_fraction_Lipschitz_exponent;
        std::vector<boost::shared_ptr<pdat::CellVariable<double> > > d_Harten_volume_fraction_Lipschitz_exponent;
        
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
        
};

#endif /* MULTIRESOLUTION_TAGGER_HPP */
