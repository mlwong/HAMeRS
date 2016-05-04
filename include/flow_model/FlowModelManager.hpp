#ifndef FLOW_MODEL_MANAGER_HPP
#define FLOW_MODEL_MANAGER_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Dimension.h"

#include "flow_model/FlowModels.hpp"

#include "equation_of_state/EquationOfStateIdealGas.hpp"
#include "flow_model/boundary_conditions/Euler/EulerBoundaryConditions.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructor.hpp"
#include "flow_model/refinement_taggers/gradient_tagger/GradientTagger.hpp"
#include "flow_model/refinement_taggers/multiresolution_tagger/MultiresolutionTagger.hpp"
#include "flow_model/initial_conditions/InitialConditions.hpp"
#include "integrator/RungeKuttaLevelIntegrator.hpp"

#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorLLF.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorFirstOrderHLLC.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWENO-JS5-LLF.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWENO-CU6-M2-LLF.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWCNS-JS5-HLLC-HLL.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWCNS-CU6-M2-HLLC-HLL.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWCNS-HW6-HLLC-HLL.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWCNS-HW6-LD-HLLC-HLL.hpp"
#include "flow_model/convective_flux_reconstructor/ConvectiveFluxReconstructorWCNS-Test.hpp"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

class FlowModelManager:
    public appu::VisDerivedDataStrategy
{
    public:
        FlowModelManager(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const std::string& flow_model_str,
            int& num_species);
        
        /*
         * Get the flow model used.
         */
        const FLOW_MODEL&
        getFlowModel() const
        {
            return d_flow_model;
        }
        
        /*
         * Get the number of equations.
         */
        const int&
        getNumberOfEquations() const
        {
            return d_num_eqn;
        }
        
        /*
         * Get the number of ghost cells.
         */
        const hier::IntVector&
        getNumberOfGhostCells() const
        {
            if (d_num_ghosts == hier::IntVector::getZero(d_dim))
            {
                TBOX_ERROR(d_object_name
                    << ": getNumberOfGhostCells()\n"
                    << "The number of ghost cells is not computed yet.\n"
                    << "The flux reconstructors may have not initialized yet."
                    << std::endl);
            }
            
            return d_num_ghosts;
        }
        
        /*
         * Initialize d_equation_of_state.
         */
        void
        initializeEquationOfState(
            const std::string& equation_of_state_string,
            const boost::shared_ptr<tbox::Database>& equation_of_state_db,
            boost::shared_ptr<EquationOfState>& equation_of_state);
        
        /*
         * Initialize d_conv_flux_reconstructor.
         * The function also computes the number of ghost cells required.
         */
        void
        initializeConvectiveFluxReconstructor(
            const std::string& shock_capturing_scheme_str,
            const boost::shared_ptr<tbox::Database>& shock_capturing_scheme_db,
            boost::shared_ptr<ConvectiveFluxReconstructor>& conv_flux_reconstructor,
            boost::shared_ptr<pdat::FaceVariable<double> >& convective_flux,
            boost::shared_ptr<pdat::CellVariable<double> >& source);
        
        /*
         * Initialize d_initial_conditions.
         */
        void
        initializeInitialConditions(
            const std::string& project_name,
            boost::shared_ptr<InitialConditions>& initial_conditions);
        
        /*
         * Initialize d_Euler_boundary_conditions.
         */
        void
        initializeEulerBoundaryConditions(
            const std::string& project_name,
            const boost::shared_ptr<tbox::Database>& Euler_boundary_conditions_db,
            const bool& Euler_boundary_conditions_db_is_from_restart,
            boost::shared_ptr<EulerBoundaryConditions>& Euler_boundary_conditions);
        
        /*
         * Initialize d_gradient_tagger.
         */
        void
        initializeGradientTagger(
            const boost::shared_ptr<tbox::Database>& gradient_tagger_db,
            boost::shared_ptr<GradientTagger>& gradient_tagger);
        
        /*
         * Initialize d_multiresolution_tagger.
         */
        void
        initializeMultiresolutionTagger(
            const boost::shared_ptr<tbox::Database>& multiresolution_tagger_db,
            boost::shared_ptr<MultiresolutionTagger>& multiresolution_tagger);
        
        /*
         * Compute and set the number of ghost cells needed.
         */
        void
        setNumberOfGhostCells();
        
        /*
         * Register the conservative variables.
         */
        void
        registerConservativeVariables(
            RungeKuttaLevelIntegrator* integrator);
        
        /*
         * Register the temporary variables used in refinement taggers.
         */
        void
        registerRefinementTaggerVariables(
            RungeKuttaLevelIntegrator* integrator);
        
        /*
         * Set the plotting context.
         */
        void
        setPlotContext(
            const boost::shared_ptr<hier::VariableContext>& plot_context);
        
        /*
         * Register the plotting quantities.
         */
#ifdef HAVE_HDF5
        void
        registerPlotQuantities(
            RungeKuttaLevelIntegrator* integrator,
            const boost::shared_ptr<appu::VisItDataWriter>& visit_writer);
#endif
        
        /*
         * Compute the stable time increment for patch using a CFL
         * condition and return the computed dt.
         */
        double
        computeStableDtOnPatch(
            hier::Patch& patch,
            const bool initial_time,
            const double dt_time,
            const boost::shared_ptr<hier::VariableContext>& data_context);
        
        /*
         * Create a vector of pointers to conservative variables with the
         * given data context.
         */
        std::vector<double*>
        createConservativeVariableVector(
            hier::Patch& patch,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const bool fill_zero = false);
        
        /*
         * Update the mass fraction/volume fraction of the last species
         * in the conservative variable vector.
         */
        void
        updateConservativeVariableVector(
            hier::Patch& patch,
            std::vector<double*>& conservative_variable_vector);
        
        /*
         * Put the d_flow_model into the restart database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Compute derived plot quantities registered with the VisIt data
         * writers from data that is maintained on each patch in the
         * hierarchy.
         */
        bool
        packDerivedDataIntoDoubleBuffer(
            double* buffer,
            const hier::Patch& patch,
            const hier::Box& region,
            const std::string& variable_name,
            int depth_id,
            double simulation_time) const;
        
        /*
         * Print all characteristics of d_flow_model.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Print all data statistics for the time-dependent variables.
         */
        void
        printDataStatistics(
            std::ostream& os,
            const boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy) const;
        
    private:
        /*
         * Initialize the number of ghost cells and boost::shared_ptr of the variables
         * in d_conv_flux_reconstructor.
         */
        void
        setVariablesForConvectiveFluxReconstructor(
            boost::shared_ptr<pdat::FaceVariable<double> >& convective_flux,
            boost::shared_ptr<pdat::CellVariable<double> >& source);
        
        /*
         * Initialize boost::shared_ptr of the variables
         * in d_initial_conditions.
         */
        void
        setVariablesForInitialConditions();
        
        /*
         * Initialize boost::shared_ptr of the variables
         * in d_Euler_boundary_conditions.
         */
        void
        setVariablesForEulerBoundaryConditions();
        
        /*
         * Initialize the number of ghost cells and boost::shared_ptr of the variables
         * in d_gradient_tagger.
         */
        void
        setVariablesForGradientTagger();
        
        /*
         * Initialize the number of ghost cells and boost::shared_ptr of the variables
         * in d_multiresolution_tagger.
         */
        void
        setVariablesForMultiresolutionTagger();
        
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
         * Number of ghost cells for time-dependent variables.
         */
        hier::IntVector d_num_ghosts;
        
        /*
         * Flow model.
         */
        FLOW_MODEL d_flow_model;
        
        /*
         * Number of equations.
         */
        int d_num_eqn;
        
        /*
         * Number of species.
         */
        const int d_num_species;
        
        /*
         * boost::shared_ptr to EquationOfState.
         */
        boost::shared_ptr<EquationOfState> d_equation_of_state;
        
        /*
         * boost::shared_ptr to the convective flux reconstructor.
         */
        boost::shared_ptr<ConvectiveFluxReconstructor> d_conv_flux_reconstructor;
        
        /*
         * boost::shared_ptr to InitialConditions.
         */
        boost::shared_ptr<InitialConditions> d_initial_conditions;
        
        /*
         * boost::shared_ptr to EulerBoundaryConditions.
         */
        boost::shared_ptr<EulerBoundaryConditions> d_Euler_boundary_conditions;
        
        /*
         * boost::shared_ptr to GradientTagger.
         */
        boost::shared_ptr<GradientTagger> d_gradient_tagger;
        
        /*
         * boost::shared_ptr to MultiresolutionTagger.
         */
        boost::shared_ptr<MultiresolutionTagger> d_multiresolution_tagger;
        
        /*
         * boost::shared_ptr to the plotting context.
         */
        boost::shared_ptr<hier::VariableContext> d_plot_context;
        
        /*
         * boost::shared_ptr to conservative variables.
         * Solution state is represented by conservative variables,
         * density, momentum, total_energy per unit volume etc.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_partial_density;
        boost::shared_ptr<pdat::CellVariable<double> > d_momentum;
        boost::shared_ptr<pdat::CellVariable<double> > d_total_energy;
        boost::shared_ptr<pdat::CellVariable<double> > d_mass_fraction;
        boost::shared_ptr<pdat::CellVariable<double> > d_volume_fraction;
        
        /*
         * Boolean to determine whether d_equation_of_state is initialized.
         */
        bool d_equation_of_state_initialized;
        
};

#endif /* FLOW_MODEL_MANAGER_HPP */
