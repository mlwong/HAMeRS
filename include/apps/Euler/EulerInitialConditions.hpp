#ifndef EULER_INITIAL_CONDITIONS_HPP
#define EULER_INITIAL_CONDITIONS_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Dimension.h"

#include "flow/flow_models/FlowModels.hpp"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

class EulerInitialConditions
{
    public:
        EulerInitialConditions(
            const std::string& object_name,
            const std::string& project_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const FLOW_MODEL::TYPE& flow_model_type,
            const int& num_species):
                d_object_name(object_name),
                d_project_name(project_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_flow_model_type(flow_model_type),
                d_num_species(num_species),
                d_variables_set(false)
        {}
        
        /*
         * Set the cell variables.
         */
        void
        setVariables(
            const std::vector<boost::shared_ptr<pdat::CellVariable<double> > >& conservative_variables)
        {
            switch (d_flow_model_type)
            {
                case FLOW_MODEL::SINGLE_SPECIES:
                {
                    d_density = conservative_variables[0];
                    d_momentum = conservative_variables[1];
                    d_total_energy = conservative_variables[2];
                    
                    d_variables_set = true;
                    
                    break;
                }
                case FLOW_MODEL::FOUR_EQN_CONSERVATIVE:
                {
                    d_partial_density = conservative_variables[0];
                    d_momentum = conservative_variables[1];
                    d_total_energy = conservative_variables[2];
                    
                    d_variables_set = true;
                    
                    break;
                }
                case FLOW_MODEL::FIVE_EQN_ALLAIRE:
                {
                    d_partial_density = conservative_variables[0];
                    d_momentum = conservative_variables[1];
                    d_total_energy = conservative_variables[2];
                    d_volume_fraction = conservative_variables[3];
                    
                    break;
                }
            }
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
         * Flow model type.
         */
        const FLOW_MODEL::TYPE d_flow_model_type;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
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

#endif /* EULER_INITIAL_CONDITIONS_HPP */
