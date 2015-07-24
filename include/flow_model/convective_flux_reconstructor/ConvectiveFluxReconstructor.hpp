#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

#include "flow_model/FlowModels.hpp"
#include "flow_model/equation_of_state/EquationOfStateIdealGas.hpp"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

class ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructor(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geom,
            const FLOW_MODEL& flow_model,
            const int& num_eqn,
            const int& num_species,
            const boost::shared_ptr<EquationOfState>& equation_of_state,
            const boost::shared_ptr<tbox::Database>& shock_capturing_scheme_db):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geom),
                d_num_conv_ghosts(hier::IntVector::getZero(d_dim)),
                d_num_ghosts(hier::IntVector::getZero(d_dim)),
                d_flow_model(flow_model),
                d_num_eqn(num_eqn),
                d_num_species(num_species),
                d_equation_of_state(equation_of_state),
                d_shock_capturing_scheme_db(shock_capturing_scheme_db),
                d_density(NULL),
                d_partial_density(NULL),
                d_momentum(NULL),
                d_total_energy(NULL),
                d_mass_fraction(NULL),
                d_volume_fraction(NULL),
                d_set_variables(false)
        {}
        
        /*
         * Get the number of ghost cells needed by the convective flux
         * reconstructor.
         */
        hier::IntVector
        getConvertiveFluxNumberOfGhostCells(void) const
        {
            return d_num_conv_ghosts;
        }
        
        /*
         * Set the cell variables if single-species flow model is chosen.
         */
        void
        setVariablesForSingleSpecies(
            const hier::IntVector& num_ghosts,
            const boost::shared_ptr<pdat::CellVariable<double> > density,
            const boost::shared_ptr<pdat::CellVariable<double> > momentum,
            const boost::shared_ptr<pdat::CellVariable<double> > total_energy,
            const boost::shared_ptr<pdat::FaceVariable<double> > convective_flux,
            const boost::shared_ptr<pdat::CellVariable<double> > source)
        {
            d_num_ghosts = num_ghosts;
            
            d_density = density;
            d_momentum = momentum;
            d_total_energy = total_energy;
            d_convective_flux = convective_flux;
            d_source = source;
            
            if (d_flow_model != SINGLE_SPECIES)
            {            
                TBOX_ERROR(d_object_name
                           << ": "
                           << "setVariablesForSingleSpecies() shouldn't be used."
                           << std::endl);
            }
            
            d_set_variables = true;
        }
        
        /*
         * Set the cell variables if four-equation multi-species flow model
         * by Shyue is chosen.
         */
        void
        setVariablesForFourEqnShyue(
            const hier::IntVector& num_ghosts,
            const boost::shared_ptr<pdat::CellVariable<double> > density,
            const boost::shared_ptr<pdat::CellVariable<double> > momentum,
            const boost::shared_ptr<pdat::CellVariable<double> > total_energy,
            const boost::shared_ptr<pdat::CellVariable<double> > mass_fraction,
            const boost::shared_ptr<pdat::FaceVariable<double> > convective_flux,
            const boost::shared_ptr<pdat::CellVariable<double> > source)
        {
            d_num_ghosts = num_ghosts;
            
            d_density = density;
            d_momentum = momentum;
            d_total_energy = total_energy;
            d_mass_fraction = mass_fraction;
            d_convective_flux = convective_flux;
            d_source = source;
            
            if (d_flow_model != FOUR_EQN_SHYUE)
            {            
                TBOX_ERROR(d_object_name
                           << ": "
                           << "setVariablesForFourEqnShyue() shouldn't be used."
                           << std::endl);
            }
            
            d_set_variables = true;
        }
        
        /*
         * Set the cell variables if five-equation multi-species flow model
         * by Allaire is chosen.
         */
        void
        setVariablesForFiveEqnAllaire(
            const hier::IntVector& num_ghosts,
            const boost::shared_ptr<pdat::CellVariable<double> > partial_density,
            const boost::shared_ptr<pdat::CellVariable<double> > momentum,
            const boost::shared_ptr<pdat::CellVariable<double> > total_energy,
            const boost::shared_ptr<pdat::CellVariable<double> > volume_fraction,
            const boost::shared_ptr<pdat::FaceVariable<double> > convective_flux,
            const boost::shared_ptr<pdat::CellVariable<double> > source)
        {
            d_num_ghosts = num_ghosts;
            
            d_partial_density = partial_density;
            d_momentum = momentum;
            d_total_energy = total_energy;
            d_volume_fraction = volume_fraction;
            d_convective_flux = convective_flux;
            d_source = source;
            
            if (d_flow_model != FIVE_EQN_ALLAIRE)
            {            
                TBOX_ERROR(d_object_name
                           << ": "
                           << "setVariablesForFiveEqnAllaire() shouldn't be used."
                           << std::endl);
            }
            
            d_set_variables = true;
        }
        
        /*
         * Print all characteristics of the convective flux reconstruction class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the convective flux reconstruction class
         * into the restart database.
         */
        virtual void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the convective fluxes and sources due to hyperbolization
         * of the equtions.
         */
        virtual void
        computeConvectiveFluxAndSource(
            hier::Patch& patch,
            const double time,
            const double dt,
            const boost::shared_ptr<hier::VariableContext> data_context) = 0;
    
    protected:
        /*
         * The object name is used for error/warning reporting and also as a
         * string label for restart database entries.
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
         * Number of ghost cells needed by the shock capturing scheme.
         */
        hier::IntVector d_num_conv_ghosts;
        
        /*
         * Number of ghost cells for time-independent variables.
         */
        hier::IntVector d_num_ghosts;
        
        /*
         * Flow model.
         */
        const FLOW_MODEL d_flow_model;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * boost::shared_ptr to EquationOfState.
         */
        const boost::shared_ptr<EquationOfState> d_equation_of_state;
        
        /*
         * boost::shared_ptr to database of the shock capturing scheme.
         */
        const boost::shared_ptr<tbox::Database> d_shock_capturing_scheme_db;
        
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
         * boost::shared_ptr to convective flux variable vector.
         */
        boost::shared_ptr<pdat::FaceVariable<double> > d_convective_flux;
        
        /*
         * boost::shared_ptr to source variable vector.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_source;
        
        /*
         * Boolean to determine where proper variables are initialized.
         */
        bool d_set_variables;
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_HPP */
