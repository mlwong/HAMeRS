#ifndef INITIAL_CONDITIONS_HPP
#define INITIAL_CONDITIONS_HPP

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

using namespace SAMRAI;

class InitialConditions
{
    public:
        InitialConditions(
            const std::string& object_name,
            const std::string& project_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const FLOW_MODEL& flow_model,
            const int& num_species,
            const boost::shared_ptr<EquationOfState>& equation_of_state):
                d_object_name(object_name),
                d_project_name(project_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_flow_model(flow_model),
                d_num_species(num_species),
                d_equation_of_state(equation_of_state),
                d_variables_set(false)
        {}
        
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
         * Set the data on the patch interior to some initial values,
         * depending on the flow problems and flow models.
         */
        void
        initializeDataOnPatch(
            hier::Patch& patch,
            const double data_time,
            const bool initial_time,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
    private:
        /*
         * The object name is used for error/warning reporting.
         */
        const std::string d_object_name;
        
        /*
         * Name of the project.
         */
        std::string d_project_name;
        
        /*
         * Problem dimension.
         */
        const tbox::Dimension d_dim;
        
        /*
         * boost::shared_ptr to the grid geometry.
         */
        const boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
        
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
};

#endif /* INITIAL_CONDITIONS_HPP */
