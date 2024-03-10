#ifndef HYPERVISCOSITY_OPERATOR_HPP
#define HYPERVISCOSITY_OPERATOR_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "flow/flow_models/FlowModels.hpp"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

#include <string>
#include <vector>

using namespace SAMRAI;

class HyperviscosityOperator
{
    public:
        HyperviscosityOperator(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const HAMERS_SHARED_PTR<FlowModel>& flow_model,
            const HAMERS_SHARED_PTR<tbox::Database>& hyperviscosity_operator_db);
        
        virtual ~HyperviscosityOperator() {}
        
        /*
         * Get the number of ghost cells needed by the hyperviscosity operator.
         */
        hier::IntVector
        getHyperviscosityOperatorNumberOfGhostCells() const
        {
            return d_num_hyperviscosity_op_ghosts;
        }
        
        /*
         * Print all characteristics of the hyperviscosity operator class.
         */
        virtual void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the hyperviscosity operator into the restart database.
         */
        virtual void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * Perform the hyperviscosity operator on a patch.
         */
        void
        performHyperviscosityOperatorOnPatch(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<hier::CoarseFineBoundary> coarse_fine_bdry,
            const HAMERS_SHARED_PTR<pdat::SideVariable<Real> >& variable_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_source,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double time,
            const double dt,
            const int RK_step_number);
        
    private:
        /*
         * Perform the hyperviscosity operator on a patch using source form.
         */
        void
        performHyperviscosityOperatorOnPatchSourceForm(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> > source,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double dt,
            const hier::Box& domain) const;
        
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
         * Number of ghost cells needed by the hyperviscosity operator.
         */
        hier::IntVector d_num_hyperviscosity_op_ghosts;
        
        /*
         * Number of equations.
         */
        const int d_num_eqn;
        
        /*
         * Flow model.
         */
        const HAMERS_SHARED_PTR<FlowModel> d_flow_model;
        
        /*
         * HAMERS_SHARED_PTR to database of the hyperviscosity operator.
         */
        const HAMERS_SHARED_PTR<tbox::Database> d_hyperviscosity_operator_db;
        
        /*
         * Scheme parameters.
         */
        int d_lap_order;
        int d_accuracy_order;
        bool d_use_flux_form;
        Real d_coeff;
        
        std::vector<Real> d_coeffs_node;
        std::vector<Real> d_coeffs_midpoint;
        Real d_prefactor;
};


#endif /* HYPERVISCOSITY_OPERATOR_HPP */
