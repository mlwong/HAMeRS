#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Serializable.h"

#include "algs/integrator/RungeKuttaLevelIntegrator.hpp"
#include "algs/patch_strategy/RungeKuttaPatchStrategy.hpp"
#include "apps/Navier-Stokes/NavierStokesBoundaryConditions.hpp"
#include "apps/Navier-Stokes/NavierStokesInitialConditions.hpp"
#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructorManager.hpp"
#include "flow/diffusive_flux_reconstructors/DiffusiveFluxReconstructorManager.hpp"
#include "flow/flow_models/FlowModelManager.hpp"
#include "flow/refinement_taggers/GradientTagger.hpp"
#include "flow/refinement_taggers/MultiresolutionTagger.hpp"
#include "flow/refinement_taggers/ValueTagger.hpp"

#include "boost/shared_ptr.hpp"
#include <string>
#include <vector>

using namespace SAMRAI;

/**
 * The NavierStokes class provides routines for a sample application code that
 * solves the Navier-Stokes equations of gas dynamics.  This code illustrates the
 * manner in which a code employing the standard Berger/Oliger AMR algorithm for
 * explicit hydrodynamics can be used in the SAMRAI framework.  This class is
 * derived from the RungeKuttaPatchStrategy abstract base class which defines the
 * bulk of the interface between the hyperbolic Runge-Kutta intergration algorithm
 * modified from SAMRAI's implementation of algs::HyperbolicPatchStrategy and the
 * numerical routines specific to NavierStokes.  In particular, this class provides
 * routines which maybe applied to any patch in an AMR patch hierarchy.
 *
 */

class NavierStokes:
    public RungeKuttaPatchStrategy,
    public tbox::Serializable
{
    public:
        /**
         * The constructor for NavierStokes sets default parameters for the
         * NavierStokes model.  Specifically, it allocates the variables that
         * represent the state of the solution.  The constructor also registers
         * this object for restart with the restart manager using the object name.
         *
         * After default values are set, this routine calls getFromRestart()
         * if execution from a restart file is specified.  Finally,
         * getFromInput() is called to read values from the given input
         * database (potentially overriding those found in the restart file).
         */
        NavierStokes(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<tbox::Database>& input_db,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry);
        
        /**
         * Destructor of NavierStokes.
         */
        ~NavierStokes();
        
        /**
         * Register NavierStokes model variables with RungeKuttaLevelIntegrator
         * according to variable registration function provided by the integrator.
         * In other words, variables are registered according to their role
         * in the integration process (e.g., time-dependent, flux, etc.).
         * This routine also registers variables for plotting with the
         * Vis writer.
         */
        void
        registerModelVariables(
            RungeKuttaLevelIntegrator* integrator);
        
        /**
         * Set up parameters in the load balancer object (owned by the gridding
         * algorithm) if needed.  The NavierStokes model allows non-uniform load
         * balancing to be used based on the input file parameter called
         * "use_nonuniform_workload".  The default case is to use uniform
         * load balancing (i.e., use_nonuniform_workload == false).  For
         * illustrative and testing purposes, when non-uniform load balancing is
         * turned on, a weight of one will be applied to every grid cell.  This
         * should produce an identical patch configuration to the uniform load
         * balance case.
         */
        void
        setupLoadBalancer(
            RungeKuttaLevelIntegrator* integrator,
            mesh::GriddingAlgorithm* gridding_algorithm);
        
        /**
         * Set the data on the patch interior to some initial values, epending
         * on the input parameters and numerical routines.  If the "initial_time"
         * flag is false, indicating that the routine is called after a regridding
         * step, the routine does nothing.
         */
        void
        initializeDataOnPatch(
            hier::Patch& patch,
            const double data_time,
            const bool initial_time);
        
        /**
         * Compute the stable time increment for patch using a CFL condition and
         * return the computed dt.
         */
        double
        computeStableDtOnPatch(
            hier::Patch& patch,
            const bool initial_time,
            const double dt_time);
        
        /**
         * Compute time integral of hyperbolic fluxes to be used in finite
         * difference for patch Runge-Kutta integration.
         * 
         * If the equations are in non-conservative form, the TIME INTEGRALS
         * of the correspoinding sources are also computed after the equation
         * is hyperolized.
         *
         * The finite difference used to update the integrated quantities
         * through the Runge-Kutta steps is implemented in the advanceSingleStep()
         * routine.
         */
        void
        computeHyperbolicFluxesAndSourcesOnPatch(
            hier::Patch& patch,
            const double time,
            const double dt,
            const int RK_step_number);
        
        /**
         * Advance a single Runge-Kutta step. Conservative differencing is
         * implemented here by using the hyperbolic fluxes and sources computed
         * in computeHyperbolicFluxesAndSourcesOnPatch()
         */
        void
        advanceSingleStep(
            hier::Patch& patch,
            const double time,
            const double dt,
            const std::vector<double>& alpha,
            const std::vector<double>& beta,
            const std::vector<double>& gamma,
            const std::vector<boost::shared_ptr<hier::VariableContext> >& intermediate_context);
        
        /**
         * Correct Navier-Stokes solution variables at coarse-fine booundaries
         * by repeating conservative differencing with corrected fluxes.
         */
        void
        synchronizeHyperbolicFluxes(
            hier::Patch& patch,
            const double time,
            const double dt);
        
        /**
         * Preprocess before tagging cells using value detector.
         */
        void
        preprocessTagValueDetectorCells(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
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
        tagValueDetectorCells(
            hier::Patch& patch,
            const double regrid_time,
            const bool initial_error,
            const int tag_indx,
            const bool uses_gradient_detector_too,
            const bool uses_multiresolution_detector_too,
            const bool uses_integral_detector_too,
            const bool uses_richardson_extrapolation_too);
        
        /**
         * Tag cells for refinement using gradient detector.
         */
        void
        tagGradientDetectorCells(
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
        preprocessTagMultiresolutionDetectorCells(
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
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
        tagMultiresolutionDetectorCells(
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
        ///  are concrete implementations of functions declared in the
        ///  RefinePatchStrategy abstract base class.
        ///
        
        /**
         * Set the data in ghost cells corresponding to physical boundary
         * conditions. Specific boundary conditions are determined by
         * information specified in input file and numerical routines.
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
        ///  are concrete implementations of functions declared in the
        ///  CoarsenPatchStrategy abstract base class.  They are trivial
        ///  because this class doesn't do any pre/postprocessCoarsen.
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
         * This routine is a concrete implementation of the function
         * declared in the tbox::Serializable abstract base class.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /**
         * Register a VisIt data writer so this class will write
         * plot files that may be postprocessed with the VisIt
         * visualization tool.
         */
#ifdef HAVE_HDF5
        void registerVisItDataWriter(const boost::shared_ptr<appu::VisItDataWriter>& viz_writer);
#endif
        
        /**
         * This routine is a concrete implementation of the virtual function
         * in the base class appu::VisDerivedDataStrategy.  It computes derived
         * plot quantities registered with the VisIt data writers from data
         * that is maintained on each patch in the hierarchy.  In particular, it
         * writes the plot quantity identified by the string variable name to
         * the specified double buffer on the patch in the given region.  The
         * depth_id integer argument indicates which entry in the "depth" of the
         * vector is being written; for a scalar quantity, this may be ignored.
         * For a vector quantity, it may be used to compute the quantity at the
         * particular depth (e.g. mom[depth_id] = rho * vel[depth_id]).  The
         * boolean return value specifies whether or not derived data exists on
         * the patch.  Generally, this will be TRUE.  If the
         * packDerivedDataIntoDoubleBuffer data does NOT exist on the patch,
         * return FALSE.
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
        ///  The following routines are specific to the NavierStokes class and
        ///  are not declared in any base class.
        ///
        
        /**
         * Print all data members for NavierStokes class.
         */
        void printClassData(std::ostream& os) const;
        
        /**
         * Print all data statistics for NavierStokes class.
         */
        void
        printDataStatistics(
            std::ostream& os,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy) const;
        
        void
        printErrorStatistics(
            std::ostream& os,
            const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy) const;

        /*
         * Set the plotting context.
         */
        void
        setPlotContext(
            const boost::shared_ptr<hier::VariableContext>& plot_context)
        {
            d_plot_context = plot_context;
        }

    private:
        /*
         * These private member functions read data from input and restart.
         * When beginning a run from a restart file, all data members are read
         * from the restart file.  If the boolean flag is true when reading
         * from input, some restart values may be overridden by those in the
         * input file.
         *
         * An assertion results if the database pointer is null.
         */
        void
        getFromInput(
            const boost::shared_ptr<tbox::Database>& input_db,
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
         * We cache pointers to the grid geometry and Vis data writers
         * to set up initial data, set physical boundary conditions,
         * and register plot variables. We also cache a pointer to the
         * plot context passed to the variable registration routine.
         */
        const boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
        
#ifdef HAVE_HDF5
        boost::shared_ptr<appu::VisItDataWriter> d_visit_writer;
#endif
        
        /*
         * Data iterms used for nonuniform load balance, if used.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_workload_variable;
        int d_workload_data_id;
        bool d_use_nonuniform_workload;
        
        /*
         * Number of ghost cells for time-independent variables.
         */
        hier::IntVector d_num_ghosts;
        
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
         * Number of species.
         */
        int d_num_species;
        
        /*
         * boost::shared_ptr to FlowModel and its database.
         */
        boost::shared_ptr<FlowModel> d_flow_model;
        boost::shared_ptr<tbox::Database> d_flow_model_db;
        
        /*
         * boost::shared_ptr to the ConvectiveFluxReconstructor and its database.
         */
        boost::shared_ptr<ConvectiveFluxReconstructor> d_convective_flux_reconstructor;
        boost::shared_ptr<tbox::Database> d_convective_flux_reconstructor_db;
        
        /*
         * boost::shared_ptr to the DiffusiveFluxReconstructor and its database.
         */
        boost::shared_ptr<DiffusiveFluxReconstructor> d_diffusive_flux_reconstructor;
        boost::shared_ptr<tbox::Database> d_diffusive_flux_reconstructor_db;
        
        /*
         * boost::shared_ptr to NavierStokesInitialConditions.
         */
        boost::shared_ptr<NavierStokesInitialConditions> d_Navier_Stokes_initial_conditions;
        
        /*
         * boost::shared_ptr to NavierStokesBoundaryConditions and its database.
         */
        boost::shared_ptr<NavierStokesBoundaryConditions> d_Navier_Stokes_boundary_conditions;
        boost::shared_ptr<tbox::Database> d_Navier_Stokes_boundary_conditions_db;
        bool d_Navier_Stokes_boundary_conditions_db_is_from_restart;
        
        /*
         * boost::shared_ptr to ValueTagger and its database.
         */
        boost::shared_ptr<ValueTagger> d_value_tagger;
        boost::shared_ptr<tbox::Database> d_value_tagger_db;
        
        /*
         * boost::shared_ptr to GradientTagger and its database.
         */
        boost::shared_ptr<GradientTagger> d_gradient_tagger;
        boost::shared_ptr<tbox::Database> d_gradient_tagger_db;
        
        /*
         * boost::shared_ptr to MultiresolutionTagger and its database.
         */
        boost::shared_ptr<MultiresolutionTagger> d_multiresolution_tagger;
        boost::shared_ptr<tbox::Database> d_multiresolution_tagger_db;
        
        /*
         * boost::shared_ptr to FlowModelManager.
         */
        boost::shared_ptr<FlowModelManager> d_flow_model_manager;
        
        /*
         * boost::shared_ptr to ConvectiveFluxReconstructorManager.
         */
        boost::shared_ptr<ConvectiveFluxReconstructorManager> d_convective_flux_reconstructor_manager;
        
        /*
         * boost::shared_ptr to DiffusiveFluxReconstructorManager.
         */
        boost::shared_ptr<DiffusiveFluxReconstructorManager> d_diffusive_flux_reconstructor_manager;
        
        /*
         * boost::shared_ptr to face variable of convective flux.
         */
        boost::shared_ptr<pdat::FaceVariable<double> > d_variable_convective_flux;
        
        /*
         * boost::shared_ptr to face variable of diffusive flux.
         */
        boost::shared_ptr<pdat::FaceVariable<double> > d_variable_diffusive_flux;
        
        /*
         * boost::shared_ptr to cell variable of source terms.
         */
        boost::shared_ptr<pdat::CellVariable<double> > d_variable_source;
        
        /*
         * boost::shared_ptr to the plotting context.
         */
        boost::shared_ptr<hier::VariableContext> d_plot_context;
        
        /*
         * Timers.
         */
        static boost::shared_ptr<tbox::Timer> t_init;
        static boost::shared_ptr<tbox::Timer> t_compute_dt;
        static boost::shared_ptr<tbox::Timer> t_compute_hyperbolicfluxes;
        static boost::shared_ptr<tbox::Timer> t_advance_steps;
        static boost::shared_ptr<tbox::Timer> t_synchronize_hyperbloicfluxes;
        static boost::shared_ptr<tbox::Timer> t_setphysbcs;
        static boost::shared_ptr<tbox::Timer> t_tagvalue;
        static boost::shared_ptr<tbox::Timer> t_taggradient;
        static boost::shared_ptr<tbox::Timer> t_tagmultiresolution;
};

#endif /* NAVIER_STOKES_HPP */
