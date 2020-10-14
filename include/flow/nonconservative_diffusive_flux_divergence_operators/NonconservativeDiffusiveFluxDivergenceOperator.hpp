#ifndef NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_HPP
#define NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_HPP

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

class NonconservativeDiffusiveFluxDivergenceOperator
{
    public:
        NonconservativeDiffusiveFluxDivergenceOperator(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& nonconservative_diffusive_flux_divergence_operator_db):
                d_object_name(object_name),
                d_dim(dim),
                d_grid_geometry(grid_geometry),
                d_num_diff_ghosts(hier::IntVector::getZero(d_dim)),
                d_num_eqn(num_eqn),
                d_flow_model(flow_model),
                d_nonconservative_diffusive_flux_divergence_operator_db(
                    nonconservative_diffusive_flux_divergence_operator_db)
        {}
        
        virtual ~NonconservativeDiffusiveFluxDivergenceOperator() {}
        
        /*
         * Get the number of ghost cells needed by the non-conservative diffusive flux divergence
         * operator.
         */
        hier::IntVector
        getNonconservativeDiffusiveFluxDivergenceOperatorNumberOfGhostCells(void) const
        {
            return d_num_diff_ghosts;
        }
        
        /*
         * Print all characteristics of the non-conservative diffusive flux divergence operator class.
         */
        virtual void
        printClassData(std::ostream& os) const = 0;
        
        /*
         * Put the characteristics of the non-conservative diffusive flux divergence operator class
         * into the restart database.
         */
        virtual void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const = 0;
        
        /*
         * Compute the non-conservative diffusive flux divergence on a patch.
         */
        virtual void
        computeNonconservativeDiffusiveFluxDivergenceOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::CellVariable<double> >& variable_diffusive_flux_divergence,
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
         * Number of ghost cells needed by the non-conservative diffusive flux divergence operator.
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
         * HAMERS_SHARED_PTR to database of the non-conservative diffusive flux divergence operatore.
         */
        const HAMERS_SHARED_PTR<tbox::Database> d_nonconservative_diffusive_flux_divergence_operator_db;
        
};

#endif /* NONCONSERVATIVE_DIFFUSIVE_FLUX_DIVERGENCE_OPERATOR_HPP */
