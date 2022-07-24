#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModels.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructor(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const FLOW_MODEL::TYPE& flow_model_type,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_conv_ghosts(hier::IntVector::getZero(d_dim)),
                d_num_eqn(num_eqn),
                d_flow_model_type(flow_model_type),
                d_flow_model(flow_model),
                d_convective_flux_reconstructor_db(convective_flux_reconstructor_db)
        {}
        
        virtual ~ConvectiveFluxReconstructor() {}
        
        /*
         * Get the number of ghost cells needed by the convective flux
         * reconstructor.
         */
        hier::IntVector
        getConvectiveFluxNumberOfGhostCells(void) const
        {
            return d_num_conv_ghosts;
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
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the convective flux and source due to splitting of convective term on a patch.
         */
        virtual void
        computeConvectiveFluxAndSourceOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::SideVariable<double> >& variable_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_source,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double time,
            const double dt,
            const int RK_step_number) = 0;
    
    protected:
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
         * Number of ghost cells needed by the convective flux reconstructor.
         */
        hier::IntVector d_num_conv_ghosts;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * Flow model type.
         */
        const FLOW_MODEL::TYPE d_flow_model_type;
        
        /*
         * Flow model.
         */
        const HAMERS_SHARED_PTR<FlowModel> d_flow_model;
        
        /*
         * HAMERS_SHARED_PTR to database of the convective_flux_reconstructor.
         */
        const HAMERS_SHARED_PTR<tbox::Database> d_convective_flux_reconstructor_db;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_HPP */
