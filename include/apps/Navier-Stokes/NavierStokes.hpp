#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include "HAMeRS_config.hpp"

#include "HAMeRS_memory.hpp"

#include "algs/integrator/RungeKuttaLevelIntegrator.hpp"
#include "algs/patch_strategy/RungeKuttaPatchStrategy.hpp"
#include "apps/Navier-Stokes/NavierStokesBoundaryConditions.hpp"
#include "apps/Navier-Stokes/NavierStokesErrorStatistics.hpp"
#include "apps/Navier-Stokes/NavierStokesInitialConditions.hpp"
#include "extn/visit_data_writer/ExtendedVisItDataWriter.hpp"
#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructorManager.hpp"
#include "flow/diffusive_flux_reconstructors/DiffusiveFluxReconstructorManager.hpp"
#include "flow/nonconservative_diffusive_flux_divergence_operators/NonconservativeDiffusiveFluxDivergenceOperatorManager.hpp"
#include "flow/flow_models/FlowModelManager.hpp"
#include "flow/refinement_taggers/GradientTagger.hpp"
#include "flow/refinement_taggers/MultiresolutionTagger.hpp"
#include "flow/refinement_taggers/ValueTagger.hpp"

#include "SAMRAI/appu/VisDerivedDataStrategy.h"
// #include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Serializable.h"

#include <string>
#include <vector>

using namespace SAMRAI;

/**
 * The NavierStokes class provides routines for a sample application code that solves the Navier-
 * Stokes equations of gas dynamics. This code illustrates the manner in which a code employing
 * the standard Berger/Oliger AMR algorithm for explicit hydrodynamics can be used in the SAMRAI
 * framework. This class is derived from the RungeKuttaPatchStrategy abstract base class which defines
 * the bulk of the interface between the hyperbolic Runge-Kutta intergration algorithm modified from
 * SAMRAI's implementation of algs::HyperbolicPatchStrategy and the numerical routines specific to
 * NavierStokes. In particular, this class provides routines which maybe applied to any patch in an
 * AMR patch hierarchy.
 *
 */

class NavierStokes:
    public RungeKuttaPatchStrategy,
    public tbox::Serializable
{
    public:
        /**
         * The constructor for Euler sets default parameters for the Euler model. Specifically, it
         * allocates the variables that represent the state of the solution. The constructor also
         * registers this object for restart with the restart manager using the object name.
         *
         * After default values are set, this routine calls getFromRestart() if execution from a
         * restart file is specified. Finally, getFromInput() is called to read values from the given
         * input database (potentially overriding those found in the restart file).
         */
        NavierStokes(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            const HAMERS_SHARED_PTR<geom::CartesianGridGeometry>& grid_geometry,
            const std::string& stat_dump_filename = "");
        
        /**
         * Destructor of NavierStokes.
         */
        ~NavierStokes();
        
        /**
         * Register NavierStokes model variables with RungeKuttaLevelIntegrator according to variable
         * registration function provided by the integrator. In other words, variables are registered
         * according to their role in the integration process (e.g., time-dependent, flux, etc.).
         * This routine also registers variables for plotting with the VisIt writer.
         */
        void
        registerModelVariables(
            RungeKuttaLevelIntegrator* integrator);
        
        /**
         * Set up parameters in the load balancer object (owned by the gridding algorithm) if needed.
         * The NavierStokes model allows non-uniform load balancing to be used based on the input file
         * parameter called "use_nonuniform_workload". The default case is to use uniform load
         * balancing (i.e., use_nonuniform_workload == false). For illustrative and testing purposes,
         * when non-uniform load balancing is turned on, a weight of one will be applied to every
         * grid cell. This should produce an identical patch configuration to the uniform load balance
         * case.
         */
        void
        setupLoadBalancer(
            RungeKuttaLevelIntegrator* integrator,
            mesh::GriddingAlgorithm* gridding_algorithm);
        
        /**
         * Set the data on the patch interior to some initial values, epending on the input parameters
         * and numerical routines. If the "initial_time" flag is false, indicating that the routine
         * is called after a regridding step, the routine does nothing.
         */
        void
        initializeDataOnPatch(
            hier::Patch& patch,
            const double data_time,
            const bool initial_time);
        
        /**
         * Get the number of spectral radiuses.
         */
        int
        getNumberOfSpectralRadiuses() const
        {
            return d_dim.getValue();
        }
        
        /**
         * Compute the spectral radiuses and local stable time increment for patch using a CFL condition
         * and return them.
         */
        std::vector<double>
        computeSpectralRadiusesAndStableDtOnPatch(
            hier::Patch& patch,
            const bool initial_time,
            const double dt_time);
        
        /**
         * Compute time integral of convective fluxes to be used in finite difference for patch Runge-
         * Kutta integration.
         * 
         * If the convective term is split into non-conservative form, the split terms are computed
         * as source terms.
         *
         * The finite difference used to update the integrated quantities through the Runge-Kutta
         * steps is implemented in the advanceSingleStepOnPatch() routine.
         */
        void
        computeFluxesAndSourcesOnPatch(
            hier::Patch& patch,
            const double time,
            const double dt,
            const int RK_step_number,
            const HAMERS_SHARED_PTR<hier::VariableContext>& data_context = HAMERS_SHARED_PTR<hier::VariableContext>());
        
        /**
         * Advance a single Runge-Kutta step. Conservative differencing is implemented here by using
         * the fluxes and sources computed in computeFluxesAndSourcesOnPatch().
         */
        void
        advanceSingleStepOnPatch(
            hier::Patch& patch,
            const double time,
            const double dt,
            const std::vector<double>& alpha,
            const std::vector<double>& beta,
            const std::vector<double>& gamma,
            const std::vector<HAMERS_SHARED_PTR<hier::VariableContext> >& intermediate_context);
        
        /**
         * Correct Navier-Stokes solution variables at coarse-fine booundaries by repeating conservative
         * differencing with corrected fluxes.
         */
        void
        synchronizeFluxes(
            hier::Patch& patch,
            const double time,
            const double dt);
        
        /**
         * Preprocess before tagging cells using value detector.
         */
        void
        preprocessTagCellsValueDetector(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const double regrid_time,
            const bool initial_error,
            const bool uses_gradient_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Tag cells for refinement using value detector.
         */
        void
        tagCellsOnPatchValueDetector(
            hier::Patch& patch,
            const double regrid_time,
            const bool initial_error,
            const int tag_indx,
            const bool uses_gradient_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Preprocess before tagging cells using gradient detector.
         */
        void
        preprocessTagCellsGradientDetector(
           const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
           const int level_number,
           const double regrid_time,
           const bool initial_error,
           const bool uses_value_detector_too,
           const bool uses_multiresolution_detector_too,
           const bool uses_integral_detector_too,
           const bool uses_richardson_extrapolation_too);
        
        /**
         * Tag cells for refinement using gradient detector.
         */
        void
        tagCellsOnPatchGradientDetector(
            hier::Patch& patch,
            const double regrid_time,
            const bool initial_error,
            const int tag_indx,
            const bool uses_value_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Preprocess before tagging cells using multiresolution detector.
         */
        void
        preprocessTagCellsMultiresolutionDetector(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const int level_number,
            const double regrid_time,
            const bool initial_error,
            const bool uses_value_detector_too,
            const bool uses_gradient_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Tag cells for refinement using multiresolution detector.
         */
        void
        tagCellsOnPatchMultiresolutionDetector(
            hier::Patch& patch,
            const double regrid_time,
            const bool initial_error,
            const int tag_indx,
            const bool uses_value_detector_too,
            const bool uses_gradient_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        //@{
        //! @name Required implementations of RungeKuttaPatchStrategy pure virtuals.
        
        ///
        ///  The following routines:
        ///
        ///      setPhysicalBoundaryConditions(),
        ///      getRefineOpStencilWidth(),
        ///      preprocessRefine(),
        ///      postprocessRefine()
        ///
        ///  are concrete implementations of functions declared in the RefinePatchStrategy abstract
        ///  base class.
        ///
        
        /**
         * Set the data in ghost cells corresponding to physical boundary conditions. Specific
         * boundary conditions are determined by information specified in input file and numerical
         * routines.
         */
        void
        setPhysicalBoundaryConditions(
            hier::Patch& patch,
            const double fill_time,
            const hier::IntVector& ghost_width_to_fill);
        
        /**
         * Return stencil width of conservative linear interpolation operations.
         */
        hier::IntVector
        getRefineOpStencilWidth(
            const tbox::Dimension& dim) const
        {
           return hier::IntVector::getOne(dim);
        }
        
        void
        preprocessRefine(
            hier::Patch& fine,
            const hier::Patch& coarse,
            const hier::Box& fine_box,
            const hier::IntVector& ratio)
        {
            NULL_USE(fine);
            NULL_USE(coarse);
            NULL_USE(fine_box);
            NULL_USE(ratio);
        }
        
        void
        postprocessRefine(
            hier::Patch& fine,
            const hier::Patch& coarse,
            const hier::Box& fine_box,
            const hier::IntVector& ratio)
        {
            NULL_USE(fine);
            NULL_USE(coarse);
            NULL_USE(fine_box);
            NULL_USE(ratio);
        }
        
        ///
        ///  The following routines:
        ///
        ///      getCoarsenOpStencilWidth(),
        ///      preprocessCoarsen()
        ///      postprocessCoarsen()
        ///
        ///  are concrete implementations of functions declared in the CoarsenPatchStrategy abstract
        ///  base class. They are trivial because this class doesn't do any pre/postprocessCoarsen.
        ///
        
        /**
         * Return stencil width of conservative averaging operations.
         */
        hier::IntVector
        getCoarsenOpStencilWidth(
            const tbox::Dimension& dim) const
        {
            return hier::IntVector::getZero(dim);
        }
        
        void
        preprocessCoarsen(
            hier::Patch& coarse,
            const hier::Patch& fine,
            const hier::Box& coarse_box,
            const hier::IntVector& ratio)
        {
            NULL_USE(coarse);
            NULL_USE(fine);
            NULL_USE(coarse_box);
            NULL_USE(ratio);
        }
        
        void
        postprocessCoarsen(
            hier::Patch& coarse,
            const hier::Patch& fine,
            const hier::Box& coarse_box,
            const hier::IntVector& ratio)
        {
            NULL_USE(coarse);
            NULL_USE(fine);
            NULL_USE(coarse_box);
            NULL_USE(ratio);
        }
        
        //@}
        
        /**
         * Write state of NavierStokes object to the given database for restart.
         *
         * This routine is a concrete implementation of the function declared in the tbox::Serializable
         * abstract base class.
         */
        void
        putToRestart(
            const HAMERS_SHARED_PTR<tbox::Database>& restart_db) const;
        
        /**
         * Register a VisIt data writer so this class will write plot files that may be postprocessed
         * with the VisIt visualization tool.
         */
#ifdef HAVE_HDF5
        void registerVisItDataWriter(const HAMERS_SHARED_PTR<ExtendedVisItDataWriter>& viz_writer);
#endif
        
        /**
         * This routine is a concrete implementation of the virtual function in the base class
         * appu::VisDerivedDataStrategy. It computes derived plot quantities registered with the
         * VisIt data writers from data that is maintained on each patch in the hierarchy. In
         * particular, it writes the plot quantity identified by the string variable name to the
         * specified double buffer on the patch in the given region. The depth_id integer argument
         * indicates which entry in the "depth" of the vector is being written; for a scalar quantity,
         * this may be ignored. For a vector quantity, it may be used to compute the quantity at the
         * particular depth (e.g. mom[depth_id] = rho * vel[depth_id]). The boolean return value
         * specifies whether or not derived data exists on the patch. Generally, this will be TRUE.
         * If the packDerivedDataIntoDoubleBuffer data does NOT exist on the patch, return FALSE.
         */
        bool
        packDerivedDataIntoDoubleBuffer(
            double* buffer,
            const hier::Patch& patch,
            const hier::Box& region,
            const std::string& variable_name,
            int depth_id,
            double simulation_time) const;
        
        ///
        ///  The following routines are specific to the NavierStokes class and are not declared in
        ///  any base class.
        ///
        
        /**
         * Print all data members for the class.
         */
        void printClassData(std::ostream& os) const;
        
        /**
         * Print data statistics (max/min conservative variables).
         */
        void
        computeAndOutputMonitoringDataStatistics(
            std::ostream& os,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const int step_num,
            const double time) const;
        
        void
        printErrorStatistics(
            std::ostream& os,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double time) const;
        
        /**
         * Compute variables for computing the statistics of data.
         */
        void
        computeStatisticsVariables(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy);
        
        /**
         * Filter variables for computing the statistics of data.
         */
        void
        filterStatisticsVariables(
            const int level,
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy);
        
        /**
         * Output the header of monitoring statistics.
         */
        void outputHeaderMonitoringStatistics();
        
        /**
         * Output the header of statistics.
         */
        void
        outputHeaderStatistics();
        
        /**
         * Compute the statistics of data.
         */
        void
        computeDataStatistics(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double statistics_data_time);
        
        /**
         * Output the statistics of data.
         */
        void
        outputDataStatistics(
            const HAMERS_SHARED_PTR<hier::PatchHierarchy>& patch_hierarchy,
            const double output_time);
        
        /**
         * Get object of storing ensemble statistics.
         */
        HAMERS_SHARED_PTR<EnsembleStatistics>
        getEnsembleStatistics();
        
        /**
         * Set object of storing ensemble statistics.
         */
        void
        setEnsembleStatistics(const HAMERS_SHARED_PTR<EnsembleStatistics> ensemble_statistics);
        
        /*
         * Set the plotting context.
         */
        void
        setPlotContext(
            const HAMERS_SHARED_PTR<hier::VariableContext>& plot_context)
        {
            d_plot_context = plot_context;
        }
        
    private:
        /*
         * These private member functions read data from input and restart. When beginning a run
         * from a restart file, all data members are read from the restart file. If the boolean flag
         * is true when reading from input, some restart values may be overridden by those in the
         * input file.
         *
         * An assertion results if the database pointer is null.
         */
        void
        getFromInput(
            const HAMERS_SHARED_PTR<tbox::Database>& input_db,
            bool is_from_restart);
        
        void getFromRestart();
        
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
         * We cache pointers to the grid geometry and Vis data writers to set up initial data, set
         * physical boundary conditions, and register plot variables. We also cache a pointer to the
         * plot context passed to the variable registration routine.
         */
        const HAMERS_SHARED_PTR<geom::CartesianGridGeometry> d_grid_geometry;
        
        /*
         * Name of file output that contains monitoring statistics of data.
         */
        const std::string d_monitoring_stat_dump_filename;
        
        /*
         * Name of file output that contains statistics of data.
         */
        const std::string d_stat_dump_filename;
        
#ifdef HAVE_HDF5
        HAMERS_SHARED_PTR<ExtendedVisItDataWriter> d_visit_writer;
#endif
        
        /*
         * Data iterms used for nonuniform load balance, if used.
         */
        HAMERS_SHARED_PTR<pdat::CellVariable<double> > d_workload_variable;
        int d_workload_data_id;
        bool d_use_nonuniform_workload;
        
        /*
         * A string variable to describe the flow model used.
         */
        std::string d_flow_model_str;
        
        /*
         * A string variable to describe the convective flux reconstructor used.
         */
        std::string d_convective_flux_reconstructor_str;
        
        /*
         * A string variable to describe the diffusive flux reconstructor used.
         */
        std::string d_diffusive_flux_reconstructor_str;
        
        /*
         * A string variable to describe the non-conservative diffusive flux divergence operator used.
         */
        std::string d_nonconservative_diffusive_flux_divergence_operator_str;
        
        /*
         * Number of species.
         */
        int d_num_species;
        
        /*
         * HAMERS_SHARED_PTR to FlowModel and its database.
         */
        HAMERS_SHARED_PTR<FlowModel> d_flow_model;
        HAMERS_SHARED_PTR<tbox::Database> d_flow_model_db;
        
        /*
         * HAMERS_SHARED_PTR to the ConvectiveFluxReconstructor and its database.
         */
        HAMERS_SHARED_PTR<ConvectiveFluxReconstructor> d_convective_flux_reconstructor;
        HAMERS_SHARED_PTR<tbox::Database> d_convective_flux_reconstructor_db;
        
        /*
         * HAMERS_SHARED_PTR to the DiffusiveFluxReconstructor and its database.
         */
        HAMERS_SHARED_PTR<DiffusiveFluxReconstructor> d_diffusive_flux_reconstructor;
        HAMERS_SHARED_PTR<tbox::Database> d_diffusive_flux_reconstructor_db;
        
        /*
         * HAMERS_SHARED_PTR to the NonconservativeDiffusiveFluxDivergenceOperator and its database.
         */
        HAMERS_SHARED_PTR<NonconservativeDiffusiveFluxDivergenceOperator>
            d_nonconservative_diffusive_flux_divergence_operator;
        HAMERS_SHARED_PTR<tbox::Database> d_nonconservative_diffusive_flux_divergence_operator_db;
        
        /*
         * Boolean to determine whether to use conservative or non-conservative form of diffusive flux.
         */
        bool d_use_conservative_form_diffusive_flux;
        
        /*
         * HAMERS_SHARED_PTR to NavierStokesInitialConditions.
         */
        HAMERS_SHARED_PTR<NavierStokesInitialConditions> d_Navier_Stokes_initial_conditions;
        
        /*
         * HAMERS_SHARED_PTR to database of initial conditions.
         */
        HAMERS_SHARED_PTR<tbox::Database> d_Navier_Stokes_initial_conditions_db;
        
        /*
         * HAMERS_SHARED_PTR to NavierStokesBoundaryConditions and its database.
         */
        HAMERS_SHARED_PTR<NavierStokesBoundaryConditions> d_Navier_Stokes_boundary_conditions;
        HAMERS_SHARED_PTR<tbox::Database> d_Navier_Stokes_boundary_conditions_db;
        bool d_Navier_Stokes_boundary_conditions_db_is_from_restart;
        
        /*
         * HAMERS_SHARED_PTR to NavierStokesErrorStatistics.
         */
        HAMERS_SHARED_PTR<NavierStokesErrorStatistics> d_Navier_Stokes_error_statistics;
        
        /*
         * HAMERS_SHARED_PTR to ValueTagger and its database.
         */
        HAMERS_SHARED_PTR<ValueTagger> d_value_tagger;
        HAMERS_SHARED_PTR<tbox::Database> d_value_tagger_db;
        
        /*
         * HAMERS_SHARED_PTR to GradientTagger and its database.
         */
        HAMERS_SHARED_PTR<GradientTagger> d_gradient_tagger;
        HAMERS_SHARED_PTR<tbox::Database> d_gradient_tagger_db;
        
        /*
         * HAMERS_SHARED_PTR to MultiresolutionTagger and its database.
         */
        HAMERS_SHARED_PTR<MultiresolutionTagger> d_multiresolution_tagger;
        HAMERS_SHARED_PTR<tbox::Database> d_multiresolution_tagger_db;
        
        /*
         * HAMERS_SHARED_PTR to FlowModelManager.
         */
        HAMERS_SHARED_PTR<FlowModelManager> d_flow_model_manager;
        
        /*
         * HAMERS_SHARED_PTR to ConvectiveFluxReconstructorManager.
         */
        HAMERS_SHARED_PTR<ConvectiveFluxReconstructorManager> d_convective_flux_reconstructor_manager;
        
        /*
         * HAMERS_SHARED_PTR to DiffusiveFluxReconstructorManager.
         */
        HAMERS_SHARED_PTR<DiffusiveFluxReconstructorManager> d_diffusive_flux_reconstructor_manager;
        
        /*
         * HAMERS_SHARED_PTR to NonconservativeDiffusiveFluxDivergenceOperatorManager.
         */
        HAMERS_SHARED_PTR<NonconservativeDiffusiveFluxDivergenceOperatorManager>
            d_nonconservative_diffusive_flux_divergence_operator_manager;
        
        /*
         * HAMERS_SHARED_PTR to side variable of convective flux.
         */
        HAMERS_SHARED_PTR<pdat::SideVariable<double> > d_variable_convective_flux;
        
        /*
         * HAMERS_SHARED_PTR to side variable of diffusive flux.
         */
        HAMERS_SHARED_PTR<pdat::SideVariable<double> > d_variable_diffusive_flux;
        
        /*
         * HAMERS_SHARED_PTR to cell variable of diffusive flux divergence.
         */
        HAMERS_SHARED_PTR<pdat::CellVariable<double> > d_variable_diffusive_flux_divergence;
        
        /*
         * HAMERS_SHARED_PTR to cell variable of source terms.
         */
        HAMERS_SHARED_PTR<pdat::CellVariable<double> > d_variable_source;
        
        /*
         * HAMERS_SHARED_PTR to the plotting context.
         */
        HAMERS_SHARED_PTR<hier::VariableContext> d_plot_context;
        
        /*
         * Timers.
         */
        static HAMERS_SHARED_PTR<tbox::Timer> t_init;
        static HAMERS_SHARED_PTR<tbox::Timer> t_compute_dt;
        static HAMERS_SHARED_PTR<tbox::Timer> t_compute_fluxes_sources;
        static HAMERS_SHARED_PTR<tbox::Timer> t_advance_step;
        static HAMERS_SHARED_PTR<tbox::Timer> t_synchronize_fluxes;
        static HAMERS_SHARED_PTR<tbox::Timer> t_setphysbcs;
        static HAMERS_SHARED_PTR<tbox::Timer> t_tagvalue;
        static HAMERS_SHARED_PTR<tbox::Timer> t_taggradient;
        static HAMERS_SHARED_PTR<tbox::Timer> t_tagmultiresolution;
        
};

#endif /* NAVIER_STOKES_HPP */
