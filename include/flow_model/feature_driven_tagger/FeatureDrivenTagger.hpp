#ifndef FEATURE_DRIVEN_TAGGER_HPP
#define FEATURE_DRIVEN_TAGGER_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Dimension.h"

#include "equation_of_state/EquationOfStateIdealGas.hpp"
#include "flow_model/FlowModels.hpp"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

class FeatureDrivenTagger
{
    public:
        FeatureDrivenTagger(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const hier::IntVector& num_ghosts,
            const FLOW_MODEL& flow_model,
            const int& num_species,
            const boost::shared_ptr<EquationOfState>& equation_of_state,
            const boost::shared_ptr<tbox::Database>& feature_driven_tagger_db);
        
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
         * Print all characteristics of the feature driven tagger class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the feature driven tagger into the restart
         * database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Tag cells for refinement using different feature-based detectors.
         */
        void
        tagCells(
            hier::Patch& patch,
            const double regrid_time,
            const bool initial_error,
            const bool uses_richardson_extrapolation_too,
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
         * Refinement criteria and tolerances for the feature-based detector.
         */
        std::vector<std::string> d_refinement_criteria;
        
        double d_density_Jameson_tol;
        double d_pressure_Jameson_tol;
        double d_shock_Ducros_tol;
        double d_shock_Larsson_tol;
        
        double d_vorticity_Q_criterion_tol;
};

#endif /* FEATURE_DRIVEN_TAGGER_HPP */
