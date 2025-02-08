#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_HPP

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

namespace WENO_INTERP
{
    enum TYPE { WENO5Z,
                WENO6LD };
}

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
            const HAMERS_SHARED_PTR<tbox::Database>& convective_flux_reconstructor_db);
        
        virtual ~ConvectiveFluxReconstructor() {}
        
        /*
         * Get the number of ghost cells needed by the convective flux
         * reconstructor.
         */
        hier::IntVector
        getConvectiveFluxNumberOfGhostCells() const
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
            const int level_number,
            const HAMERS_SHARED_PTR<hier::CoarseFineBoundary>& coarse_fine_bdry,
            const HAMERS_SHARED_PTR<pdat::SideVariable<Real> >& variable_convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellVariable<Real> >& variable_source,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const double time,
            const double dt,
            const int RK_step_number) = 0;
    
    protected:
        /*
         * Put the characteristics of the base convective flux reconstruction class into the restart database.
         */
        void
        putToRestartBase(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /*
         * (Old) Compute the convective flux and source due to splitting using shock-capturing scheme.
         */
        virtual void
        computeConvectiveFluxAndSourceOnPatchShockCapturingOld(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& source_scratch,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const hier::Box& domain,
            const double dt,
            const bool use_shock_capturing,
            const bool use_interface_capturing) const;
        
        /*
         * Compute the convective flux and source due to splitting using shock-capturing scheme.
         */
        virtual void
        computeConvectiveFluxAndSourceOnPatchShockCapturing(
            hier::Patch& patch,
            const HAMERS_SHARED_PTR<pdat::SideData<Real> >& convective_flux,
            const HAMERS_SHARED_PTR<pdat::CellData<Real> >& source_scratch,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context,
            const hier::Box& domain,
            const double dt,
            const bool use_shock_capturing,
            const bool use_interface_capturing) const;
        
        /*
         * Perform WENO interpolation.
         */
        void
        performWENOInterpolation(
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& variables_minus,
            std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > >& variables_plus,
            const std::vector<std::vector<HAMERS_SHARED_PTR<pdat::SideData<Real> > > >& variables,
            const std::vector<hier::Box>& domains) const;
        
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
         * HAMERS_SHARED_PTR to database of the convective flux reconstructor.
         */
        const HAMERS_SHARED_PTR<tbox::Database> d_convective_flux_reconstructor_db;
        
        /*
         * Forms of equations.
         */
        std::vector<EQN_FORM::TYPE> d_eqn_form;
        bool d_has_advective_eqn_form;
        
        /*
         * Constants for shock- and interface-capturing scheme.
         */
        
        Real d_threshold_sensor_shock;
        Real d_threshold_sensor_interface;
        
        const int d_num_ghosts_shock_interface_capturing;
        
        WENO_INTERP::TYPE d_weno_interp;
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_HPP */
