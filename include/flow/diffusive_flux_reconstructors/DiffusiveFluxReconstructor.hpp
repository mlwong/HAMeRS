#ifndef DIFFUSIVE_FLUX_RECONSTRUCTOR_HPP
#define DIFFUSIVE_FLUX_RECONSTRUCTOR_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModels.hpp"
#include "flow/flow_models/FlowModelSubgridScaleModel.hpp"

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

class DiffusiveFluxReconstructor
{
    public:
        DiffusiveFluxReconstructor(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& diffusive_flux_reconstructor_db):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_diff_ghosts(hier::IntVector::getZero(d_dim)),
                d_num_eqn(num_eqn),
                d_flow_model(flow_model),
                d_diffusive_flux_reconstructor_db(diffusive_flux_reconstructor_db)
        {}
        
        virtual ~DiffusiveFluxReconstructor() {}
        
        /*
         * Get the number of ghost cells needed by the diffusive flux
         * reconstructor.
         */
        hier::IntVector
        getDiffusiveFluxNumberOfGhostCells(void) const
        {
            return d_num_diff_ghosts;
        }
        
        /*
         * Print all characteristics of the diffusive flux reconstruction class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the diffusive flux reconstruction class
         * into the restart database.
         */
        virtual void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the diffusive flux on a patch.
         */
        virtual void
        computeDiffusiveFluxOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::SideVariable<Real> >& variable_diffusive_flux,
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
         * Number of ghost cells needed by the diffusive flux reconstructor.
         */
        hier::IntVector d_num_diff_ghosts;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * Flow model.
         */
        const HAMERS_SHARED_PTR<FlowModel> d_flow_model;
        
        /*
         * HAMERS_SHARED_PTR to database of the diffusive flux reconstructor.
         */
        const HAMERS_SHARED_PTR<tbox::Database> d_diffusive_flux_reconstructor_db;
        
};

#endif /* DIFFUSIVE_FLUX_RECONSTRUCTOR_HPP */
